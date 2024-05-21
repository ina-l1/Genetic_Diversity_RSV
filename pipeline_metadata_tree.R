#Pipeline: Add important metadata (Season, Clade), Data Sorting and Tree Visualization 

library(dplyr)
library(MMWRweek)
library(ggplot2)
library(scales)
library(ape)
library(treeio)
library(ggtree)

# Sequence information for NCBI WGS sequences starting from 2014/2015 season 
rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/NCBI_rsvA_wgs_germany_2015.csv")
rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/NCBI_rsvB_wgs_germany_2015.csv")

refA <- read.csv("~/RSV/git/RSV Genetic Diversity/Refseq_rsvA_NC001803.csv") #Metadata for reference strains
refB <- data.frame(Accession = "AY353550.1")

rsvA_nextclade <- read.csv("~/RSV/git/RSV Genetic Diversity/nextclade_rsvA.csv") #Nextclade output
rsvB_nextclade <- read.csv("~/RSV/git/RSV Genetic Diversity/nextclade_rsvB.csv")

# Reformat Dates
rsvA$Collection_Date <- as.Date(rsvA$Collection_Date, format = "%m/%d/%Y")
rsvB$Collection_Date <- as.Date(rsvB$Collection_Date, format = "%m/%d/%Y")

rsvA$Collection_Year <- as.numeric(format(rsvA$Collection_Date,"%Y"))
rsvA$Collection_Month <- as.numeric(format(rsvA$Collection_Date,"%m"))

rsvB$Collection_Year <- as.numeric(format(rsvB$Collection_Date, "%Y"))
rsvB$Collection_Month <- as.numeric(format(rsvB$Collection_Date, "%m"))

# RSV epidemic season in Germany from Week 40 to Week 39 of next year (October)
rsvAMMWR <- MMWRweek(rsvA$Collection_Date)
rsvBMMWR <- MMWRweek(rsvB$Collection_Date)
rsvA$MMWRyear <- rsvAMMWR$MMWRyear
rsvA$MMWRweek <- rsvAMMWR$MMWRweek
rsvA$MMWRday <- rsvAMMWR$MMWRday
rsvB$MMWRyear <- rsvBMMWR$MMWRyear
rsvB$MMWRweek <- rsvBMMWR$MMWRweek
rsvB$MMWRday <- rsvBMMWR$MMWRday

rsvA$Collection_Season <- ifelse(rsvA$MMWRweek >= 1 & rsvA$MMWRweek <= 39, paste(rsvA$MMWRyear-1,"/",rsvA$MMWRyear, sep = ""), paste(rsvA$MMWRyear,"/",rsvA$MMWRyear+1, sep = ""))
rsvB$Collection_Season <- ifelse(rsvB$MMWRweek >= 1 & rsvB$MMWRweek <= 39, paste(rsvB$MMWRyear-1,"/",rsvB$MMWRyear, sep = ""), paste(rsvB$MMWRyear,"/",rsvB$MMWRyear+1, sep = ""))

# A and B in one table
rsvA$Type <- "A"
rsvB$Type <- "B"
rsvAB <- rbind(rsvA, rsvB)

# rsvA + RefSeqA, rsvB + RefSeqB
rsvA_ref <- full_join(rsvA, refA, by = join_by(Accession, Organism_Name, Organization, Isolate,
                                               Molecule_type, Length, Nuc_Completeness, Country, Collection_Date, GenBank_Title))
rsvB_ref <- full_join(rsvB, refB, by = join_by(Accession))

# Sort dataframes
rsvA <- rsvA[order(rsvA$Accession),]
rsvB <- rsvB[order(rsvB$Accession),]
rsvAB <- rsvAB[order(rsvAB$Accession),]
rsvA_ref <- rsvA_ref[order(rsvA_ref$Accession),]
rsvB_ref <- rsvB_ref[order(rsvB_ref$Accession),]

# Add nextclade clade information to metadata
rsvA_nextclade <- arrange(rsvA_nextclade, seqName)
rsvB_nextclade <- arrange(rsvB_nextclade, seqName)
rsvA_clades <- rsvA_nextclade[c("clade", "G_clade")]
rsvB_clades <- rsvB_nextclade[c("clade", "G_clade")]
rsvA_ref <- cbind(rsvA_ref, rsvA_clades)
rsvB_ref <- cbind(rsvB_ref, rsvB_clades)

# Export metadata
write.csv(rsvA_ref, file = "~/RSV/git/RSV Genetic Diversity/rsvA_ref_metadata.csv", row.names = FALSE)
write.csv(rsvB_ref, file = "~/RSV/git/RSV Genetic Diversity/rsvB_ref_metadata.csv", row.names = FALSE)

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

# Tree visualization: RSV A
tree_A <- read.tree("~/RSV/git/RSV Genetic Diversity/Treefiles/tree_rsvA_MAFFT_alignment.treefile") #read.nexus for NEXUS format, read.tree for Newick format

treelabel_A <- tree_A$tip.label
treelabel_A <- sort(treelabel_A)

df_label_A <- data.frame(treelabel = treelabel_A,
                       accession = rsvA_ref$Accession,
                       type = rsvA_ref$Type,
                       collection_date = rsvA_ref$Collection_Date,
                       season = rsvA_ref$Collection_Season) 

df_label_A$tip.label <- paste(df_label_A$accession, df_label_A$type, df_label_A$season, sep = "_")
df_label_A$tip.label[2] <- "NC001803.1_A_RefSeq" #Reference genome

tree_A <- rename_taxa(tree_A, data = df_label_A, key = 1, value = tip.label) #rename tip label of tree into shorter version

treeA_plot <- ggtree(tree_A)
treeA_plot

df_label_A <- relocate(df_label_A, tip.label) #move tip.label to first column for plotting

treeA_plot_new <- treeA_plot %<+% df_label_A + 
  geom_tippoint(aes(color = season), size=20) +
  geom_tiplab(
    geom = "text",
    size = 10
  ) +
  guides(color = guide_legend(title = "Season")) +
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )
treeA_plot_new

#ggsave(filename = "~/RSV/git/RSV Genetic Diversity/Treefiles/tree_rsvA_annotated.pdf", width = 300, height = 300, units = "cm", limitsize = FALSE)

# Tree visualization: RSV B

tree_B <- read.tree("~/RSV/git/RSV Genetic Diversity/Treefiles/tree_rsvB_MAFFT_alignment.treefile")

treelabel_B <- tree_B$tip.label
treelabel_B <- sort(treelabel_B)

df_label_B <- data.frame(treelabel = treelabel_B,
                         accession = rsvB_ref$Accession,
                         type = rsvB_ref$Type,
                         collection_date = rsvB_ref$Collection_Date,
                         season = rsvB_ref$Collection_Season) 

df_label_B$tip.label <- paste(df_label_B$accession, df_label_B$type, df_label_B$season, sep = "_")
df_label_B$tip.label[1] <- "AY353550.1_B_RefSeq"

tree_B <- rename_taxa(tree_B, data = df_label_B, key = 1, value = tip.label)

treeB_plot <- ggtree(tree_B)
treeB_plot

df_label_B <- relocate(df_label_B, tip.label)

treeB_plot_new <- treeB_plot %<+% df_label_B +
  geom_tippoint(aes(color = season), size = 20) +
  geom_tiplab(
    geom = "text",
    size = 10
  ) +
  guides(color = guide_legend(title = "Season")) +
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )
treeB_plot_new
