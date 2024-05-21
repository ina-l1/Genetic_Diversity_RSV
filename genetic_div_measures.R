# Calculate Genetic Diverstiy 

library(dplyr)
library(ape)
library(bio3d)
library(ggplot2)
library(plotly)
library(adegenet)
library(usedist)

# Read sequence alignment and metadata
aln_rsvA <- read.FASTA("~/RSV/git/RSV Genetic Diversity/Sequences/rsvA_MAFFT_alignment.fasta")
aln_rsvB <- read.FASTA("~/RSV/git/RSV Genetic Diversity/Sequences/rsvB_MAFFT_alignment.fasta")

aln_shannon_rsvA <- read.fasta("~/RSV/git/RSV Genetic Diversity/Sequences/rsvA_MAFFT_alignment.fasta")

# Important metadata
meta_rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/rsvA_ref_metadata.csv")
meta_rsvA[2, "Collection_Season"] <- "Ref"
meta_rsvA_short <- subset(meta_rsvA, select = c("Accession", "Type", "Collection_Date", "Collection_Season", "Collection_Month"))
meta_rsvA_short$plotlabel <- paste(meta_rsvA_short$Accession, meta_rsvA_short$Type, meta_rsvA_short$Collection_Season, meta_rsvA_short$Collection_Month, sep = "_")

meta_rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/rsvB_ref_metadata.csv")
meta_rsvB[1, "Collection_Season"] <- "Ref"
meta_rsvB_short <- subset(meta_rsvB, select = c("Accession", "Type", "Collection_Date", "Collection_Season", "Collection_Month"))
meta_rsvB_short$plotlabel <- paste(meta_rsvB_short$Accession, meta_rsvB_short$Type, meta_rsvB_short$Collection_Season, meta_rsvB_short$Collection_Month, sep = "_")

# Pairwise Distances from DNA Sequences
# Maybe use dist() for euclidian distance
dist_rsvA <- dist.dna(aln_rsvA, model = "K80", #evolutionary model 
                      variance = FALSE, #compute variances of distances
                      gamma = FALSE, #correction of distances
                      pairwise.deletion = FALSE, #delete sites with missing data
                      as.matrix = TRUE) #return results as matrix or object of class dist

dist_rsvB <- dist.dna(aln_rsvB, model = "K80", 
                      variance = FALSE, 
                      gamma = FALSE,
                      pairwise.deletion = FALSE,
                      as.matrix = TRUE)

# SNP distance with snp-dists
snp_dist_rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/snpdist_rsvA.csv")
rownames(snp_dist_rsvA) <- snp_dist_rsvA$snp.dists.0.7.0
snp_dist_rsvA <- select(snp_dist_rsvA, -snp.dists.0.7.0)

snp_dist_rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/snpdist_rsvB.csv")
rownames(snp_dist_rsvB) <- snp_dist_rsvB$snp.dists.0.7.0
snp_dist_rsvB <- select(snp_dist_rsvB, -snp.dists.0.7.0)

# Shannon Entropy
'shannon_rsvA <- entropy(aln_shannon_rsvA)
shannon_H_rsvA <- shannon_rsvA$H
shannon_H_rsvA_df <- as.data.frame(shannon_H_rsvA)
shannon_H_rsvA_df$position <- as.numeric(rownames(shannon_H_rsvA_df))
colnames(shannon_H_rsvA_df) <- c("shannon_entropy", "position")
shannon_H_rsvA_df <- relocate(shannon_H_rsvA_df, position)

shannon_plot <- ggplot(shannon_H_rsvA_df, aes(x = position, y = shannon_entropy)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 15000, 1000))
shannon_plot'

# Classic multidimensional scaling (MDS)

# Pairwise distance
mds_rsvA <- cmdscale(dist_rsvA, k = 2) #k dim
mds_rsvB <- cmdscale(dist_rsvB, k = 2)

# SNP distance
mds_snp_rsvA <- cmdscale(snp_dist_rsvA, k = 2)
mds_snp_rsvB <- cmdscale(snp_dist_rsvB, k = 2)

# Plot MDS

#Pairwise distance 
x_mds_rsvA <- mds_rsvA[, 1]
y_mds_rsvA <- mds_rsvA[, 2]

plot_mds_rsvA <- plot(x_mds_rsvA, y_mds_rsvA)

x_mds_rsvB <- mds_rsvB[, 1]
y_mds_rsvB <- mds_rsvB[, 2]

plot_mds_rsvB <- plot(x_mds_rsvB, y_mds_rsvB)

# SNP distance
x_mds_snp_rsvA <- mds_snp_rsvA[, 1]
y_mds_snp_rsvA <- mds_snp_rsvA[, 2]

plot_mds_snp_rsvA <- plot(x_mds_snp_rsvA, y_mds_snp_rsvA)

x_mds_snp_rsvB <- mds_snp_rsvB[, 1]

y_mds_snp_rsvB <- mds_snp_rsvB[, 2]

plot_mds_snp_rsvB <- plot(x_mds_snp_rsvB, y_mds_snp_rsvB)

# Plot MDS with labels 

# PAIRWISE: Create df with labels
df_mds_rsvA <- data.frame(x_axis = x_mds_rsvA, y_axis = y_mds_rsvA)
df_mds_rsvA$ID <- rownames(df_mds_rsvA)
df_mds_rsvA <- relocate(df_mds_rsvA, ID)
rownames(df_mds_rsvA) <- c(1:nrow(df_mds_rsvA))
df_mds_rsvA <- arrange(df_mds_rsvA, ID)

mds_meta_rsvA <- cbind(meta_rsvA_short, df_mds_rsvA)

df_mds_rsvB <- data.frame(x_axis = x_mds_rsvB, y_axis = y_mds_rsvB)
df_mds_rsvB$ID <- rownames(df_mds_rsvB)
df_mds_rsvB <- relocate(df_mds_rsvB, ID)
rownames(df_mds_rsvB) <- c(1:nrow(df_mds_rsvB))
df_mds_rsvB <- arrange(df_mds_rsvB, ID)

mds_meta_rsvB <- cbind(meta_rsvB_short, df_mds_rsvB)

# SNP: Create df with labels
df_mds_snp_rsvA <- data.frame(x_axis = x_mds_snp_rsvA, y_axis = y_mds_snp_rsvA)
df_mds_snp_rsvA$ID <- rownames(df_mds_snp_rsvA)
df_mds_snp_rsvA <- relocate(df_mds_snp_rsvA, ID)
rownames(df_mds_snp_rsvA) <- c(1:nrow(df_mds_snp_rsvA))
df_mds_snp_rsvA <- arrange(df_mds_snp_rsvA, ID)

df_mds_snp_rsvB <- data.frame(x_axis = x_mds_snp_rsvB, y_axis = y_mds_snp_rsvB)
df_mds_snp_rsvB$ID <- rownames(df_mds_snp_rsvB)
df_mds_snp_rsvB <- relocate(df_mds_rsvB, ID)
rownames(df_mds_snp_rsvB) <- c(1:nrow(df_mds_snp_rsvB))
df_mds_snp_rsvB <- arrange(df_mds_snp_rsvB, ID)

mds_snp_meta_rsvA <- cbind(meta_rsvA_short, df_mds_snp_rsvA)

mds_snp_meta_rsvB <- cbind(meta_rsvB_short, df_mds_snp_rsvB)

# MDS Plots

# Pairwise distance
plot_rsvA_mds_label <- ggplot(mds_meta_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) + #color = Collection_Season
  geom_point() +
  geom_text(
    label = paste(mds_meta_rsvA$Collection_Season, mds_meta_rsvA$Collection_Month, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvA_mds_label

plot_rsvB_mds_label <- ggplot(mds_meta_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_meta_rsvB$Collection_Season, mds_meta_rsvB$Collection_Month, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvB_mds_label

# SNP distance
plot_rsvA_mds_snp_label <- ggplot(mds_snp_meta_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_meta_rsvA$Collection_Season, mds_snp_meta_rsvA$Collection_Month, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvA_mds_snp_label

plot_rsvB_mds_snp_label <- ggplot(mds_snp_meta_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_meta_rsvB$Collection_Season, mds_snp_meta_rsvB$Collection_Month, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvB_mds_snp_label

# Remove Ref

# Pairwise distance
mds_noref_rsvA <- filter(mds_meta_rsvA, Collection_Season != "Ref")

plot_rsvA_mds_label_noref <- ggplot(mds_noref_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_noref_rsvA$Collection_Season, mds_noref_rsvA$Collection_Month, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS for RSV-A sequences") +
  theme(
    axis.title = element_blank()
  )
plot_rsvA_mds_label_noref

mds_noref_rsvB <- filter(mds_meta_rsvB, Collection_Season != "Ref")

plot_rsvB_mds_label_noref <- ggplot(mds_noref_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_noref_rsvB$Collection_Season, mds_noref_rsvB$Collection_Month, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS for RSV-B sequences") +
  theme(
    axis.title = element_blank()
  )
plot_rsvB_mds_label_noref

# SNP distance 
mds_snp_noref_rsvA <- filter(mds_snp_meta_rsvA, Collection_Season != "Ref")

plot_rsvA_mds_snp_label_noref <- ggplot(mds_snp_noref_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_noref_rsvA$Collection_Season, mds_snp_noref_rsvA$Collection_Month, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS SNP Distance RSV-A") +
  theme(
    axis.title = element_blank()
  )
plot_rsvA_mds_snp_label_noref

mds_snp_noref_rsvB <- filter(mds_snp_meta_rsvB, Collection_Season != "Ref")

plot_rsvB_mds_snp_label_noref <- ggplot(mds_snp_noref_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_noref_rsvB$Collection_Season, mds_snp_noref_rsvB$Collection_Month, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS SNP Distance RSV-B") +
  theme(
    axis.title = element_blank()
  )
plot_rsvB_mds_snp_label_noref

# Violin Plot 

# Pairwise distance
dist_rsvA <- dist_rsvA[order(rownames(snp_dist_rsvA)),]
distgroup_pair_rsvA <- dist_groups(dist_rsvA, meta_rsvA_short$Collection_Season)

distgroup_within_pair_rsvA <- subset(distgroup_pair_rsvA, Group1 == Group2)
distgroup_within_pair_rsvA$Season <- distgroup_within_pair_rsvA$Group1

distgroup_between_pair_rsvA <- subset(distgroup_pair_rsvA, Group1 != Group2)
distgroup_between_pair_rsvA$Between_Season <- paste(distgroup_between_pair_rsvA$Group1, 
                                                    "&",
                                                    distgroup_between_pair_rsvA$Group2,
                                                    sep = "")

plot_violin_pair_within_rsvA <- ggplot(distgroup_within_pair_rsvA, aes(x = Season, y = Distance)) +
  geom_violin()
plot_violin_pair_within_rsvA

plot_violin_pair_between_rsvA <- ggplot(distgroup_between_pair_rsvA, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_pair_between_rsvA

dist_rsvB <- dist_rsvB[order(rownames(snp_dist_rsvB)),]
distgroup_pair_rsvB <- dist_groups(dist_rsvB, meta_rsvB_short$Collection_Season)

distgroup_within_pair_rsvB <- subset(distgroup_pair_rsvB, Group1 == Group2)
distgroup_within_pair_rsvB$Within_Season <- distgroup_within_pair_rsvB$Group1

distgroup_between_pair_rsvB <- subset(distgroup_pair_rsvB, Group1 != Group2)

distgroup_between_pair_rsvB <- subset(distgroup_pair_rsvB, Group1 != Group2)
distgroup_between_pair_rsvB$Between_Season <- paste(distgroup_between_pair_rsvB$Group1, 
                                                    "&",
                                                    distgroup_between_pair_rsvB$Group2,
                                                    sep = "")

plot_violin_pair_rsvB <- ggplot(distgroup_within_pair_rsvB, aes(x = Within_Season, y = Distance)) +
  geom_violin()
plot_violin_pair_rsvB

plot_violin_pair_between_rsvB <- ggplot(distgroup_between_pair_rsvB, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_pair_between_rsvB

# SNP distance
snp_dist_rsvA <- snp_dist_rsvA[order(rownames(snp_dist_rsvA)),] #sort alphabetically
distgroup_snp_rsvA <- dist_groups(snp_dist_rsvA, meta_rsvA_short$Collection_Season)

distgroup_within_snp_rsvA <- subset(distgroup_snp_rsvA, Group1 == Group2)
distgroup_within_snp_rsvA$Season <- distgroup_within_snp_rsvA$Group1

distgroup_between_snp_rsvA <- subset(distgroup_snp_rsvA, Group1 != Group2)
distgroup_between_snp_rsvA$Between_Season <- paste(distgroup_between_snp_rsvA$Group1, 
                                                    "&",
                                                   distgroup_between_snp_rsvA$Group2,
                                                    sep = "")

plot_violin_snp_rsvA <- ggplot(distgroup_within_snp_rsvA, aes(x = Season, y = Distance)) +
  geom_violin()
plot_violin_snp_rsvA

plot_violin_snp_between_rsvA <- ggplot(distgroup_between_snp_rsvA, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_snp_between_rsvA

snp_dist_rsvB <- snp_dist_rsvB[order(rownames(snp_dist_rsvB)),] #sort alphabetically
distgroup_snp_rsvB <- dist_groups(snp_dist_rsvB, meta_rsvB_short$Collection_Season)

distgroup_within_snp_rsvB <- subset(distgroup_snp_rsvB, Group1 == Group2)
distgroup_within_snp_rsvB$Within_Season <- distgroup_within_snp_rsvB$Group1

distgroup_between_snp_rsvB <- subset(distgroup_snp_rsvB, Group1 != Group2)
distgroup_between_snp_rsvB$Between_Season <- paste(distgroup_between_snp_rsvB$Group1, 
                                                   "&",
                                                   distgroup_between_snp_rsvB$Group2,
                                                   sep = "")

plot_violin_snp_rsvB <- ggplot(distgroup_within_snp_rsvB, aes(x = Within_Season, y = Distance)) +
  geom_violin()
plot_violin_snp_rsvB

plot_violin_snp_between_rsvB <- ggplot(distgroup_between_snp_rsvB, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_snp_between_rsvB