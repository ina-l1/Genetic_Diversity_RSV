#Add important metadata (Season, Clade), Data Sorting and Tree Visualization 

library(dplyr)
library(MMWRweek)
library(ggplot2)
library(scales)
library(ape)
library(treeio)
library(ggtree)

# Metadata
RefSeq_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/RefSeq/RefSeq_rsvA.csv")
RefSeq_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/RefSeq/RefSeq_rsvB.csv")

rsvA_ref <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv")
rsvB_ref <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

rsvA_ref$Collection_Date <- as.Date(rsvA_ref$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
rsvB_ref$Collection_Date <- as.Date(rsvB_ref$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y

rsvA_ref <- subset(rsvA_ref, Accession != "NC_038235.1") ###
rsvB_ref <- subset(rsvB_ref, Accession != "NC_001781.1") ###
rsvA_ref$'EU/GER' <- ifelse(rsvA_ref$Country == "Germany", "GER", "EU") ###
rsvB_ref$'EU/GER' <- ifelse(rsvB_ref$Country == "Germany", "GER", "EU") ###

# Tree visualization: RSV A

#tree_A <- read.tree("~/Yale_Projects/Genetic_Diversity_RSV/Europe/Treefiles/tree_rsvA_MAFFT_alignment_EU.treefile") #read.nexus for NEXUS format, read.tree for Newick format
tree_A <- read.tree("~/Yale_Projects/Genetic_Diversity_RSV/BEAST/tree_rsvA_EU.tre") 

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
  geom_tippoint(aes(color = season), size=20) + #color = season
  geom_tiplab(
    geom = "text",
    size = 10
  ) +
#  scale_color_manual(values = c("GER" = "red", "EU" =  "blue")) +
  guides(color = guide_legend(title = "Season")) + #title = "Season"
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )
treeA_plot_new

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/BEAST/tree_rsvA_EU_season.pdf", width = 300, height = 1000, units = "cm", limitsize = FALSE)

# Tree visualization: RSV B

#tree_B <- read.tree("~/Yale_Projects/Genetic_Diversity_RSV/Europe/Treefiles/tree_rsvB_MAFFT_alignment_EU.treefile")
tree_B <- read.tree("~/Yale_Projects/Genetic_Diversity_RSV/BEAST/tree_rsvB_EU.tree") 

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

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/BEAST/tree_rsvB_EU.pdf", width = 300, height = 1000, units = "cm", limitsize = FALSE)
