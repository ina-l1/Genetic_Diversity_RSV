#Add important metadata (Season, Clade), Data Sorting and Tree Visualization 

library(dplyr)
library(MMWRweek)
library(ggplot2)
library(ape)
library(treeio)
library(ggtree)
library(lubridate)

##################################################

# Base directory

base_dir <- "~/Yale_Projects/Genetic_Diversity_RSV/"

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

# Tree visualization: RSV A

path_tree_A <- file.path(base_dir, "BEAST", "tree_files", "tree_rsvA_EU.trees")
tree_A <- read.tree(path_tree_A) #read.nexus for NEXUS format, read.tree for Newick format

treelabel_A <- tree_A$tip.label
treelabel_A <- sort(treelabel_A)

df_label_A <- data.frame(treelabel = treelabel_A,
                       accession = rsvA_ref$Accession,
                       type = rsvA_ref$Type,
                       collection_date = rsvA_ref$Collection_Date,
                       season = rsvA_ref$Collection_Season,
                       country = rsvA_ref$Country,
                       EUGER = rsvA_ref$`EU/GER`)

df_label_A$tip.label <- paste(df_label_A$accession, df_label_A$type, df_label_A$season, df_label_A$country, sep = "_")
#df_label_A$tip.label[which(df_label_A$accession == RefSeq_rsvA$Accession)] <- paste(RefSeq_rsvA$Accession, "A", "Ref", sep = "_") #Reference genome

tree_A <- rename_taxa(tree_A, data = df_label_A, key = 1, value = tip.label) #rename tip label of tree into shorter version

treeA_plot <- ggtree(tree_A)
treeA_plot

df_label_A <- relocate(df_label_A, tip.label) #move tip.label to first column for plotting

treeA_plot_new <- treeA_plot %<+% df_label_A + 
  geom_tippoint(aes(color = country), size=20) + #color = season
  geom_tiplab(
    geom = "text",
    size = 10
  ) +
#  scale_color_manual(values = c("GER" = "red", "EU" =  "blue")) +
  guides(color = guide_legend(title = "Country")) + #title = "Season"
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )
treeA_plot_new

path_treeA_plot_new <- file.path(base_dir, "BEAST", "tree_rsvA_EU_country.pdf")
#ggsave(filename = path_treeA_plot_new, width = 500, height = 1000, units = "cm", limitsize = FALSE)

# Tree visualization: RSV B

path_tree_B <- file.path(base_dir, "BEAST", "tree_files", "tree_rsvB_EU.trees")
tree_B <- read.tree(path_tree_B) #read.nexus for NEXUS format, read.tree for Newick format

treelabel_B <- tree_B$tip.label
treelabel_B <- sort(treelabel_B)

df_label_B <- data.frame(treelabel = treelabel_B,
                         accession = rsvB_ref$Accession,
                         type = rsvB_ref$Type,
                         collection_date = rsvB_ref$Collection_Date,
                         season = rsvB_ref$Collection_Season,
                         country = rsvB_ref$Country,
                         EUGER = rsvB_ref$`EU/GER`) 

df_label_B$tip.label <- paste(df_label_B$accession, df_label_B$type, df_label_B$season, df_label_B$country, sep = "_")
#df_label_B$tip.label[which(df_label_B$accession == RefSeq_rsvB$Accession)] <- paste(RefSeq_rsvB$Accession, "B", "Ref", sep = "_") #Reference genome

tree_B <- rename_taxa(tree_B, data = df_label_B, key = 1, value = tip.label)

treeB_plot <- ggtree(tree_B)
treeB_plot

df_label_B <- relocate(df_label_B, tip.label)

treeB_plot_new <- treeB_plot %<+% df_label_B +
  geom_tippoint(aes(color = EUGER), size = 20) +
  geom_tiplab(
    geom = "text",
    size = 10
  ) +
  scale_color_manual(values = c("GER" = "red", "EU" =  "blue")) +
  guides(color = guide_legend(title = "EU/GER")) +
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  ) 
treeB_plot_new

path_treeB_plot_new <- file.path(base_dir, "BEAST", "tree_rsvB_EU_GER.pdf")
#ggsave(filename = path_treeA_plot_new, width = 500, height = 1000, units = "cm", limitsize = FALSE)