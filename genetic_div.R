# Calculate Genetic Diversity
# Visualize using MDS
# Violin plots show genetic diversity within and between season 

library(dplyr)
library(stringr)
library(ape)
library(bio3d)
library(ggplot2)
library(plotly)
library(adegenet)
library(usedist)

# Read sequence alignment and metadata

meta_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv")
meta_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

aln_rsvA <- read.FASTA("~/Yale_Projects/Genetic_Diversity_RSV/Europe/Sequences/rsvA_MAFFT_alignment_EU.fasta")
aln_rsvB <- read.FASTA("~/Yale_Projects/Genetic_Diversity_RSV/Europe/Sequences/rsvB_MAFFT_alignment_EU.fasta")

RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

#aln_shannon_rsvA <- read.fasta("~/Yale_Projects/Genetic_Diversity_RSV/Sequences/rsvA_MAFFT_alignment.fasta")
 
# Filter important metadata
meta_rsvA$'EU/GER' <- ifelse(meta_rsvA$Country == "Germany", "GER", "EU")
meta_rsvA[which(meta_rsvA$Accession == RefSeq_rsvA), c("Collection_Season", "EU/GER", "Country")] <- "Ref"

meta_rsvA_short <- subset(meta_rsvA, select = c("Accession", "Type", "Collection_Date", "Collection_Season", "MMWRyear", "MMWRweek", "Country", "EU/GER"))
meta_rsvA_short$plotlabel <- paste(meta_rsvA_short$Accession, meta_rsvA_short$Collection_Season, meta_rsvA_short$MMWRweek, meta_rsvA_short$Country, sep = "_")

meta_rsvB$'EU/GER' <- ifelse(meta_rsvB$Country == "Germany", "GER", "EU")
meta_rsvB[which(meta_rsvB$Accession == RefSeq_rsvB), c("Collection_Season", "EU/GER", "Country")] <- "Ref"

meta_rsvB_short <- subset(meta_rsvB, select = c("Accession", "Type", "Collection_Date", "Collection_Season", "MMWRyear", "MMWRweek", "Country", "EU/GER"))
meta_rsvB_short$plotlabel <- paste(meta_rsvB_short$Accession, meta_rsvB_short$Collection_Season, meta_rsvB_short$MMWRweek, meta_rsvB_short$Country, sep = "_")

# Pairwise Distances from DNA Sequences

# Maybe use dist() for euclidian distance
dist_rsvA <- dist.dna(aln_rsvA, model = "K81", #evolutionary model 
                      variance = FALSE, #compute variances of distances
                      gamma = FALSE, #correction of distances
                      pairwise.deletion = FALSE, #delete sites with missing data
                      as.matrix = TRUE) #return results as matrix or object of class dist

dist_rsvB <- dist.dna(aln_rsvB, model = "K81", 
                      variance = FALSE, 
                      gamma = FALSE,
                      pairwise.deletion = FALSE,
                      as.matrix = TRUE)

#write.csv(dist_rsvA, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/pairdist_rsvA_EU.csv")
#write.csv(dist_rsvB, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/pairdist_rsvB_EU.csv")

# SNP distance with snp-dists

snp_dist_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU.csv")
rownames(snp_dist_rsvA) <- snp_dist_rsvA[,1]
snp_dist_rsvA <- snp_dist_rsvA[, -1]

snp_dist_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU.csv")
rownames(snp_dist_rsvB) <- snp_dist_rsvB[, 1]
snp_dist_rsvB <- snp_dist_rsvB[, -1]

# Hamming Distance

ham_dist_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_protein.csv")
rownames(ham_dist_rsvA) <- ham_dist_rsvA[, 1]
ham_dist_rsvA <- ham_dist_rsvA[, -1]

ham_dist_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_protein.csv")
rownames(ham_dist_rsvB) <- ham_dist_rsvB[, 1]
ham_dist_rsvB <- ham_dist_rsvB[, -1]

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
mds_pair_rsvA_eig <- cmdscale(dist_rsvA, k = 2, eig = TRUE) #k dim
mds_rsvA <- as.data.frame(mds_pair_rsvA_eig[1])

mds_pair_rsvB_eig <- cmdscale(dist_rsvB, k = 2, eig = TRUE)
mds_rsvB <- as.data.frame(mds_pair_rsvB_eig[1])

# SNP distance
mds_snp_rsvA_eig <- cmdscale(snp_dist_rsvA, k = 2, eig = TRUE)
mds_snp_rsvA <- as.data.frame(mds_snp_rsvA_eig[1])

mds_snp_rsvB_eig <- cmdscale(snp_dist_rsvB, k = 2, eig = TRUE)
mds_snp_rsvB <- as.data.frame(mds_snp_rsvB_eig[1])

# Hamming distance
mds_ham_rsvA_eig <- cmdscale(ham_dist_rsvA, k = 2, eig = TRUE)
mds_ham_rsvA <- as.data.frame(mds_ham_rsvA_eig[1])

mds_ham_rsvB_eig <- cmdscale(ham_dist_rsvB, k = 2, eig = TRUE)
mds_ham_rsvB <- as.data.frame(mds_ham_rsvB_eig[1])

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

# Hamming distance
x_mds_ham_rsvA <- mds_ham_rsvA[, 1]
y_mds_ham_rsvA <- mds_ham_rsvA[, 2]

plot_mds_ham_rsvA <- plot(x_mds_ham_rsvA, y_mds_ham_rsvA)

x_mds_ham_rsvB <- mds_ham_rsvB[, 1]
y_mds_ham_rsvB <- mds_ham_rsvB[, 2]

plot_mds_ham_rsvB <- plot(x_mds_ham_rsvB, y_mds_ham_rsvB)

# Calculate proportion of variance explained through MDS

## Pairwise distance
eigenvalues_pair <- mds_pair_rsvA_eig$eig
total_variance_pair <- sum(eigenvalues_pair[eigenvalues_pair > 0])
variance_explained_dim1_pair <- eigenvalues_pair[1] / total_variance_pair
variance_explained_dim2_pair <- eigenvalues_pair[2] / total_variance_pair
# Combined variance explained by the two dimensions
combined_variance_explained_pair <- variance_explained_dim1_pair + variance_explained_dim2_pair
combined_variance_explained_pair

eigenvalues_pair <- mds_pair_rsvB_eig$eig
total_variance_pair <- sum(eigenvalues_pair[eigenvalues_pair > 0])
variance_explained_dim1_pair <- eigenvalues_pair[1] / total_variance_pair
variance_explained_dim2_pair <- eigenvalues_pair[2] / total_variance_pair
# Combined variance explained by the two dimensions
combined_variance_explained_pair <- variance_explained_dim1_pair + variance_explained_dim2_pair
combined_variance_explained_pair

## SNP distance
eigenvalues_snp <- mds_snp_rsvA_eig$eig
total_variance_snp <- sum(eigenvalues_snp[eigenvalues_snp > 0])
variance_explained_dim1_snp <- eigenvalues_snp[1] / total_variance_snp
variance_explained_dim2_snp <- eigenvalues_snp[2] / total_variance_snp
# Combined variance explained by the two dimensions
combined_variance_explained_snp <- variance_explained_dim1_snp + variance_explained_dim2_snp
combined_variance_explained_snp

eigenvalues_snp <- mds_snp_rsvB_eig$eig
total_variance_snp <- sum(eigenvalues_snp[eigenvalues_snp > 0])
variance_explained_dim1_snp <- eigenvalues_snp[1] / total_variance_snp
variance_explained_dim2_snp <- eigenvalues_snp[2] / total_variance_snp
# Combined variance explained by the two dimensions
combined_variance_explained_snp <- variance_explained_dim1_snp + variance_explained_dim2_snp
combined_variance_explained_snp

## Hamming distance
eigenvalues_ham <- mds_ham_rsvA_eig$eig
total_variance_ham <- sum(eigenvalues_ham[eigenvalues_ham > 0])
variance_explained_dim1_ham <- eigenvalues_ham[1] / total_variance_ham
variance_explained_dim2_ham <- eigenvalues_ham[2] / total_variance_ham
# Combined variance explained by the two dimensions
combined_variance_explained_ham <- variance_explained_dim1_ham + variance_explained_dim2_ham

# Plot MDS with labels 

# PAIRWISE: Create df with labels
df_mds_rsvA <- data.frame(x_axis = x_mds_rsvA, y_axis = y_mds_rsvA)
df_mds_rsvA$ID <- as.numeric(rownames(df_mds_rsvA))
df_mds_rsvA <- relocate(df_mds_rsvA, ID)
rownames(df_mds_rsvA) <- c(1:nrow(df_mds_rsvA))
df_mds_rsvA <- arrange(df_mds_rsvA, ID)

mds_meta_rsvA <- cbind(meta_rsvA_short, df_mds_rsvA)

df_mds_rsvB <- data.frame(x_axis = x_mds_rsvB, y_axis = y_mds_rsvB)
df_mds_rsvB$ID <- as.numeric(rownames(df_mds_rsvB))
df_mds_rsvB <- relocate(df_mds_rsvB, ID)
rownames(df_mds_rsvB) <- c(1:nrow(df_mds_rsvB))
df_mds_rsvB <- arrange(df_mds_rsvB, ID)

mds_meta_rsvB <- cbind(meta_rsvB_short, df_mds_rsvB)

# SNP: Create df with labels
df_mds_snp_rsvA <- data.frame(x_axis = x_mds_snp_rsvA, y_axis = y_mds_snp_rsvA)
df_mds_snp_rsvA$ID <- as.numeric(rownames(df_mds_snp_rsvA))
df_mds_snp_rsvA <- relocate(df_mds_snp_rsvA, ID)
rownames(df_mds_snp_rsvA) <- c(1:nrow(df_mds_snp_rsvA))
df_mds_snp_rsvA <- arrange(df_mds_snp_rsvA, ID)

mds_snp_meta_rsvA <- cbind(meta_rsvA_short, df_mds_snp_rsvA)

df_mds_snp_rsvB <- data.frame(x_axis = x_mds_snp_rsvB, y_axis = y_mds_snp_rsvB)
df_mds_snp_rsvB$ID <- as.numeric(rownames(df_mds_snp_rsvB))
df_mds_snp_rsvB <- relocate(df_mds_snp_rsvB, ID)
rownames(df_mds_snp_rsvB) <- c(1:nrow(df_mds_snp_rsvB))
df_mds_snp_rsvB <- arrange(df_mds_snp_rsvB, ID)

mds_snp_meta_rsvB <- cbind(meta_rsvB_short, df_mds_snp_rsvB)

# Hamming: Create df with labels
df_mds_ham_rsvA <- data.frame(x_axis = x_mds_ham_rsvA, y_axis = y_mds_ham_rsvA)
df_mds_ham_rsvA$ID <- as.numeric(rownames(df_mds_ham_rsvA))
df_mds_ham_rsvA <- relocate(df_mds_ham_rsvA, ID)
rownames(df_mds_ham_rsvA) <- c(1:nrow(df_mds_ham_rsvA))
df_mds_ham_rsvA <- arrange(df_mds_ham_rsvA, ID)

mds_ham_meta_rsvA <- cbind(meta_rsvA_short, df_mds_ham_rsvA)

df_mds_ham_rsvB <- data.frame(x_axis = x_mds_ham_rsvB, y_axis = y_mds_ham_rsvB)
df_mds_ham_rsvB$ID <- as.numeric(rownames(df_mds_ham_rsvB))
df_mds_ham_rsvB <- relocate(df_mds_ham_rsvB, ID)
rownames(df_mds_ham_rsvB) <- c(1:nrow(df_mds_ham_rsvB))
df_mds_ham_rsvB <- arrange(df_mds_ham_rsvB, ID)

mds_ham_meta_rsvB <- cbind(meta_rsvB_short, df_mds_ham_rsvB)

# MDS Plots

# Pairwise distance

plot_rsvB_mds_label <- ggplot(mds_meta_rsvB, aes(x = x_axis, y = y_axis, color = Country)) + #color = Country
  geom_point() +
  theme_minimal()
plot_rsvB_mds_label

plot_rsvB_mds_label <- ggplot(mds_meta_rsvB, aes(x = x_axis, y = y_axis, color = `EU/GER`)) + #color = Collection_Season
  geom_point() +
  scale_color_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) +
  scale_fill_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) + 
  theme_minimal()
plot_rsvB_mds_label

plot_rsvB_mds_label <- ggplot(mds_meta_rsvB, aes(x = x_axis, y = y_axis, color = `EU/GER`)) + #color = Collection_Season
  geom_point() +
  geom_text(
    label = mds_meta_rsvB$plotlabel,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) +
  scale_fill_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) +
  theme_minimal()
plot_rsvB_mds_label

# SNP distance
plot_rsvB_mds_snp_label <- ggplot(mds_snp_meta_rsvB, aes(x = x_axis, y = y_axis, color = Country)) +
  geom_point() +
  theme_minimal()
plot_rsvB_mds_snp_label

plot_rsvB_mds_snp_label <- ggplot(mds_snp_meta_rsvB, aes(x = x_axis, y = y_axis, color = `EU/GER`)) +
  geom_point() +
  scale_color_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) +
  scale_fill_manual(values = c("GER" = "red", "EU" =  "blue", "Ref" = "green")) + 
  theme_minimal()
plot_rsvB_mds_snp_label

plot_rsvB_mds_snp_label <- ggplot(mds_snp_meta_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_meta_rsvB$Collection_Season, mds_snp_meta_rsvB$MMWRweek, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvB_mds_snp_label

# Hamming distance
plot_rsvA_mds_ham_label <- ggplot(mds_ham_meta_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_ham_meta_rsvA$Collection_Season, mds_ham_meta_rsvA$MMWRweek, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvA_mds_ham_label

plot_rsvB_mds_ham_label <- ggplot(mds_ham_meta_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_ham_meta_rsvB$Collection_Season, mds_ham_meta_rsvB$MMWRweek, sep = "_"),
    check_overlap = TRUE
  )
plot_rsvB_mds_ham_label

# Remove Ref

# Pairwise distance
mds_noref_rsvA <- filter(mds_meta_rsvA, Collection_Season != "Ref")

plot_rsvA_mds_label_noref <- ggplot(mds_noref_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_noref_rsvA$Collection_Season, mds_noref_rsvA$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS Pairwise Distance RSV-A") +
  theme(
    axis.title = element_blank()
  )
plot_rsvA_mds_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_pair_rsvA_noref_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

mds_noref_rsvB <- filter(mds_meta_rsvB, Collection_Season != "Ref")

plot_rsvB_mds_label_noref <- ggplot(mds_noref_rsvB, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_noref_rsvB$Collection_Season, mds_noref_rsvB$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS Pairwise Distance RSV-B") +
  theme(
    axis.title = element_blank()
  )
plot_rsvB_mds_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_pair_rsvB_noref_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

# SNP distance 
mds_snp_noref_rsvA <- filter(mds_snp_meta_rsvA, Collection_Season != "Ref")

plot_rsvA_mds_snp_label_noref <- ggplot(mds_snp_noref_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_noref_rsvA$Collection_Season, mds_snp_noref_rsvA$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS SNP Distance RSV-A") +
  theme(
    axis.title = element_blank()
  )
plot_rsvA_mds_snp_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_snp_rsvA_noref_Country_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

mds_snp_noref_rsvB <- filter(mds_snp_meta_rsvB, Collection_Season != "Ref")

plot_rsvB_mds_snp_label_noref <- ggplot(mds_snp_noref_rsvB, aes(x = x_axis, y = y_axis, color = Country)) +
  geom_point() +
  geom_text(
    label = paste(mds_snp_noref_rsvB$Collection_Season, mds_snp_noref_rsvB$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS SNP Distance RSV-B") +
  theme(
    axis.title = element_blank()
  )
plot_rsvB_mds_snp_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_snp_rsvB_noref_Country_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

# Hamming distance
mds_ham_noref_rsvA <- filter(mds_ham_meta_rsvA, Collection_Season != "Ref")

plot_rsvA_mds_ham_label_noref <- ggplot(mds_ham_noref_rsvA, aes(x = x_axis, y = y_axis, color = Collection_Season)) +
  geom_point() +
  geom_text(
    label = paste(mds_ham_noref_rsvA$Collection_Season, mds_ham_noref_rsvA$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS ham Distance RSV-A") +
  theme(
    axis.title = element_blank()
  )
plot_rsvA_mds_ham_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_ham_rsvA_noref_Country_EU.png", width = 45, height = 40, units = "cm", limitsize = FALSE)

mds_ham_noref_rsvB <- filter(mds_ham_meta_rsvB, Collection_Season != "Ref")

plot_rsvB_mds_ham_label_noref <- ggplot(mds_ham_noref_rsvB, aes(x = x_axis, y = y_axis, color = Country)) +
  geom_point() +
  geom_text(
    label = paste(mds_ham_noref_rsvB$Collection_Season, mds_ham_noref_rsvB$MMWRweek, sep = "_"),
    check_overlap = TRUE
  ) +
  labs(title = "MDS ham Distance RSV-B") +
  theme(
    axis.title = element_blank()
  )
plot_rsvB_mds_ham_label_noref

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/mds_snp_rsvB_noref_Country_EU.png", width = 45, height = 40, units = "cm", limitsize = FALSE)

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
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_pair_within_rsvA

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvA_pair_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

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
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_pair_rsvB

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvB_pair_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

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
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_snp_rsvA

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvA_snp_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

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
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_snp_rsvB

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvB_snp_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

plot_violin_snp_between_rsvB <- ggplot(distgroup_between_snp_rsvB, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_snp_between_rsvB

# Hamming distance
ham_dist_rsvA <- ham_dist_rsvA[order(rownames(ham_dist_rsvA)),] #sort alphabetically
distgroup_ham_rsvA <- dist_groups(ham_dist_rsvA, meta_rsvA_short$Collection_Season)

distgroup_within_ham_rsvA <- subset(distgroup_ham_rsvA, Group1 == Group2)
distgroup_within_ham_rsvA$Season <- distgroup_within_ham_rsvA$Group1

distgroup_between_ham_rsvA <- subset(distgroup_ham_rsvA, Group1 != Group2)
distgroup_between_ham_rsvA$Between_Season <- paste(distgroup_between_ham_rsvA$Group1, 
                                                   "&",
                                                   distgroup_between_ham_rsvA$Group2,
                                                   sep = "")

plot_violin_ham_rsvA <- ggplot(distgroup_within_ham_rsvA, aes(x = Season, y = Distance)) +
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_ham_rsvA

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvA_ham_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

plot_violin_ham_between_rsvA <- ggplot(distgroup_between_ham_rsvA, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_ham_between_rsvA


ham_dist_rsvB <- ham_dist_rsvB[order(rownames(ham_dist_rsvB)),] #sort alphabetically
distgroup_ham_rsvB <- dist_groups(ham_dist_rsvB, meta_rsvB_short$Collection_Season)

distgroup_within_ham_rsvB <- subset(distgroup_ham_rsvB, Group1 == Group2)
distgroup_within_ham_rsvB$Within_Season <- distgroup_within_ham_rsvB$Group1

distgroup_between_ham_rsvB <- subset(distgroup_ham_rsvB, Group1 != Group2)
distgroup_between_ham_rsvB$Between_Season <- paste(distgroup_between_ham_rsvB$Group1, 
                                                   "&",
                                                   distgroup_between_ham_rsvB$Group2,
                                                   sep = "")

plot_violin_ham_rsvB <- ggplot(distgroup_within_ham_rsvB, aes(x = Within_Season, y = Distance)) +
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot",
               geom = "crossbar", width = 0.05, color = "red")
plot_violin_ham_rsvB

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/violin_rsvB_ham_within_EU.png", width = 25, height = 20, units = "cm", limitsize = FALSE)

plot_violin_ham_between_rsvB <- ggplot(distgroup_between_ham_rsvB, aes(x = Between_Season, y = Distance)) +
  geom_violin()
plot_violin_ham_between_rsvB
