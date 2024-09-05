# Investigating temporal dynamics of genetic diversity in RSV subgroups during an epidemic cycle in Germany, 2014-2023

## Overview

This is the GitHub repository for the following ongoing project: <br>
*Investigating temporal dynamics of genetic diversity in RSV subgroups during an epidemic cycle in Germany, 2014-2023*

This study uses RSV genomic data to analyze the temporal dynamics of RSV within the epidemic season. For further information please refer to the abstract below.  

### Files and folder structure

**datasorting.R** Read GenBank output and create metadata tables.
**timeseries.R** Visualization of sequence and case numbers.
**genetic_div.R** Calculating distance matrices and visualize using MDS (multidimensional scaling).
**sliding_window.R** Implementation of 8-week sliding window to visualize the dynamics of genetic diversity over time.
**ref_dist.R** Pairwise distance of sequences to reference sequence in 8-week sliding window.
**phylo_tree.R** Reading and annotating time-scaled phylogenetic trees (BEAST output).
**skygrid.R** Plotting skygrid data (BEAST output).

### Installation 

All available code is written in R (v4.4.1) using Rstudio (v2024.04.2+764). 

To run code offline please adjust the file paths accordingly. 

#### Required R packages

- tidyr (1.3.1): Data manipulation
- ggplot2 (3.5.1): Plots and data visualization
- patchwork (1.2.0): Plot manipulation
- ggtree (3.12.0): Phylogenetic tree visualization and annotation
- treeio (1.28.0): Read phylogenetic trees
- lubridate (1.9.3): Handling dates 
- MMWRweek (0.1.3): Converting dates to MMWR format
- ape (5.8): Reading and analyzing genomic data
- usedist (0.4.0): Calculate distance matrix

`# Required packages
libs <- c("tidyr", "ggplot2", "patchwork", "ggtree", "treeio", "lubridate", "MMWRweek", "ape", "usedist")

# Install and load missing packages
sapply(libs, function(pkg) {
  if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg)
  library(pkg, character.only = TRUE)
})`

### Other software

Bayesian phylogenetic analysis was conducted using BEAST (Suchard MA *et al.*, 2018).

>Suchard MA, Lemey P, Baele G, Ayres DL, Drummond AJ & Rambaut A (2018) Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10 Virus Evolution 4, vey016. DOI: 10.1093/ve/vey01

SNP distance was calculated using snp-dists (v0.7.0).

>Seemann, T. (2018). Source code for snp-dists software (0.6.2). Zenodo. https://doi.org/10.5281/zenodo.1411986

Sequence data was visualized and edited using SeaView (v5.0.5) and UGENE (v50.0)

>Gouy, M. Guindon, S. & Gascuel., O. (2010) SeaView version 4 : a multiplatform graphical user interface for sequence alignment and phylogenetic tree building. Molecular Biology and Evolution 27(2):221-224. 
Galtier, N., Gouy, M. & Gautier, C. (1996) SEAVIEW and PHYLO_WIN: two graphic tools for sequence alignment and molecular phylogeny. Comput. Appl. Biosci., 12:543-548.

>Okonechnikov K, Golosova O, Fursov M, the UGENE team.  Unipro UGENE: a unified bioinformatics toolkit . Bioinformatics  2012 28: 1166-1167. doi:10.1093/bioinformatics/bts091

## Abstract

(*currently under review*)

### Background and aims

The genome of Respiratory Syncytial Virus (RSV) contains 10 genes encoding 11 proteins. RSV is divided into two antigenic subgroups, A and B, primarily based on differences in the attachment (G) glycoprotein. Both subtypes co-circulate, but one typically dominates during an epidemic season. While the genetic variability and evolutionary patterns of RSV at the genotype level have been studied, gaps remain in understanding the diversity of circulating lineages within a single season, and the mechanisms behind shifts in subgroup dominance between seasons. This study evaluates genetic diversity and population dynamics across multiple seasons to determine the spread of RSV-A and -B lineages within an epidemic cycle, with a focus on sequences from Germany.

### Methods

We analyzed 168 RSV-A and 135 RSV-B whole genome sequences (WGS) from Germany, alongside 563 RSV-A and 564 RSV-B sequences from other European countries retrieved from NCBI GenBank. Only sequences with collection dates from the 2014/15 season to the 2022/23 season were considered. Non-phylogenetic methods assessed pairwise genetic diversity using single nucleotide polymorphism (SNP) and evolutionary distances with an 8-week sliding window. Bayesian phylogenetic analysis was conducted, and time-stamped phylogenetic trees and Bayesian SkyGrids were used to assess temporal and population dynamics of RSV-A and -B. 

### Results 

For German RSV sequences, we observed a decrease in average pairwise distance of the dominant subtype (e.g. RSV-A) towards the end of the season if a shift in type dominance occurred in the following season (e.g. RSV-A to RSV-B). This pattern was also observed in sequences from Spain. RSV-B generally showed lower variation in genetic diversity than RSV-A. Phylogenetic analysis indicated that German RSV-A and RSV-B sequences from a single season typically belong to multiple clades. However, we observed that the majority of German post-COVID-19 pandemic RSV-A and RSV-B sequences each derived from a single monophyletic clade. 

### Conclusion

Examining the temporal pattern of genetic diversity within an epidemic cycle could demonstrate the circulation of RSV lineages throughout a season and reveal potential indicators for subtype dominance changes before the following season. The increased resolution achieved through this analysis provides further insight for future RSV surveillance and public health strategies. 

### Author list

Ina Li, RWTH Aachen University, Germany / Department of Epidemiology of Microbial Diseases, Yale School of Public Health, USA.
Jiye Kwon, Department of Epidemiology of Microbial Diseases, Yale School of Public Health, USA.
Verity Hill, Department of Epidemiology of Microbial Diseases, Yale School of Public Health, USA.
C. Brandon Ogbunu, Department of Ecology and Evolutionary Biology, Yale School of Medicine, USA.
Anna Matuszyńska, Computational Life Science, Department of Biology, RWTH Aachen University, Germany.
Virginia E. Pitzer, Department of Epidemiology of Microbial Diseases, Yale School of Public Health, USA.
Daniel M. Weinberger, Department of Epidemiology of Microbial Diseases, Yale School of Public Health, USA.