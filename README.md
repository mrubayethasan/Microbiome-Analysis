# Microbiome Analysis of bronchoalveolar lavage (BAL) specimens

## Overview
This repository contains the R scripts used to analyze BAL microbiome data from immunocompromised patients with pneumonia compared to healthy volunteers.

- **Analysis environment**: R (RStudio v2026.1.1.403)  
- **Key packages**:
  - vegan
  - MaAsLin2

### Analytical steps:
1. Data import and preprocessing  
2. Alpha diversity analysis (Wilcoxon Rank-Sum test)  
3. Beta diversity analysis (Bray–Curtis, Jaccard, PERMANOVA)  
4. Differential abundance testing (MaAsLin2)  
5. Visualization and figure generation  

# ==================================
# Microbiome Analysis of BAL Samples
# ==================================

# ---------------------------
# Load Libraries
# ---------------------------
library(tidyverse)
library(readxl)
library(vegan)
library(ape)
library(ggpubr)
library(reshape2)
library(Maaslin2)

# ---------------------------
# Load Data
# ---------------------------
# Replace with your own paths or use relative paths
abund <- read_excel("data/abundance_table.xlsx")
meta  <- read_excel("data/metadata.xlsx")

meta$Group <- factor(meta$Group)

# ---------------------------
# Preprocessing
# ---------------------------
tax_col <- colnames(abund)[1]

abund_mat <- abund %>%
  column_to_rownames(tax_col) %>%
  as.matrix()

# Match sample order
abund_mat <- abund_mat[, meta$SampleID]

# Relative abundance
abund_rel <- sweep(abund_mat, 2, colSums(abund_mat), "/")
abund_rel[is.na(abund_rel)] <- 0

# ---------------------------
# Taxonomy Parsing
# ---------------------------
tax_split <- str_split_fixed(rownames(abund_mat), ";", 8)

taxonomy <- data.frame(
  Phylum  = tax_split[,3],
  Genus   = tax_split[,7],
  Species = tax_split[,8]
)
rownames(taxonomy) <- rownames(abund_mat)

# =====================
# 1. Relative Abundance
# =====================

# -------- Phylum level --------
phylum_abund <- aggregate(abund_rel, by=list(taxonomy$Phylum), sum)
rownames(phylum_abund) <- phylum_abund$Group.1
phylum_abund$Group.1 <- NULL

phylum_long <- melt(phylum_abund, id.vars="row.names")
colnames(phylum_long) <- c("Phylum","SampleID","Abundance")

phylum_long <- left_join(phylum_long, meta, by="SampleID") %>%
  filter(Phylum != "Unknown")

ggplot(phylum_long, aes(SampleID, Abundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ Group, scales="free_x", space="free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

# -------- Top Genera --------
genus_abund <- aggregate(abund_rel, by=list(taxonomy$Genus), sum)
rownames(genus_abund) <- genus_abund$Group.1
genus_abund$Group.1 <- NULL

genus_long <- melt(genus_abund, id.vars="row.names")
colnames(genus_long) <- c("Genus","SampleID","Abundance")

genus_long <- left_join(genus_long, meta, by="SampleID")

top_genera <- genus_long %>%
  group_by(Genus) %>%
  summarise(Total=sum(Abundance)) %>%
  slice_max(Total, n=10) %>%
  pull(Genus)

genus_long <- genus_long %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"))

ggplot(genus_long, aes(SampleID, Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ Group, scales="free_x", space="free") +
  theme_classic()

# -------- Top Species --------
species_abund <- aggregate(abund_rel, by=list(taxonomy$Species), sum)
rownames(species_abund) <- species_abund$Group.1
species_abund$Group.1 <- NULL

species_long <- melt(species_abund, id.vars="row.names")
colnames(species_long) <- c("Species","SampleID","Abundance")

species_long <- left_join(species_long, meta, by="SampleID")

top_species <- species_long %>%
  group_by(Species) %>%
  summarise(Total=sum(Abundance)) %>%
  slice_max(Total, n=10) %>%
  pull(Species)

species_long <- species_long %>%
  mutate(Species = ifelse(Species %in% top_species, Species, "Other"))

ggplot(species_long, aes(SampleID, Abundance, fill=Species)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ Group, scales="free_x", space="free") +
  theme_classic()

# ==================
# 2. Alpha Diversity
# ==================

alpha_df <- data.frame(
  SampleID = colnames(abund_mat),
  Observed = specnumber(t(abund_mat)),
  Shannon  = diversity(t(abund_mat), "shannon"),
  Simpson  = diversity(t(abund_mat), "simpson"),
  Chao1    = estimateR(t(abund_mat))["S.chao1", ]
)

alpha_df <- left_join(alpha_df, meta, by="SampleID")

alpha_long <- pivot_longer(alpha_df,
                           cols = c(Observed, Shannon, Simpson, Chao1),
                           names_to="Metric",
                           values_to="Value")

ggplot(alpha_long, aes(Group, Value, fill=Group)) +
  geom_boxplot() +
  geom_jitter(width=0.15) +
  facet_wrap(~Metric, scales="free_y") +
  stat_compare_means(method="wilcox.test") +
  theme_classic()

# =================
# 3. Beta Diversity
# =================

# Bray-Curtis
dist_bray <- vegdist(t(abund_rel), method="bray")
pcoa_bray <- pcoa(dist_bray)

# Jaccard
dist_jaccard <- vegdist(t(abund_rel > 0), method="jaccard")
pcoa_jaccard <- pcoa(dist_jaccard)

# PERMANOVA
adonis2(dist_bray ~ Group, data=meta)
adonis2(dist_jaccard ~ Group, data=meta)

# =====================================
# 4. Differential Abundance (MaAsLin2)
# =====================================

# Prepare species-level data
abund$Species <- sapply(strsplit(abund[[1]], ";"), tail, 1)

abund_species <- abund %>%
  select(-1) %>%
  mutate(Species = abund$Species) %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Species")

abund_rel_species <- sweep(abund_species, 2, colSums(abund_species), "/")

meta_da <- meta %>%
  column_to_rownames("SampleID")

fit <- Maaslin2(
  input_data = t(abund_rel_species),
  input_metadata = meta_da,
  output = "results/maaslin2",
  fixed_effects = "Group",
  normalization = "NONE",
  transform = "LOG"
)

# ==========================
# 5. Top 20 Species Heatmap
# ==========================

df_long <- abund %>%
  mutate(Species = sub(".*;", "", abund[[1]])) %>%
  pivot_longer(-c(1, Species), names_to="SampleID", values_to="Abundance")

top20 <- df_long %>%
  group_by(Species) %>%
  summarise(Total=sum(Abundance)) %>%
  slice_max(Total, n=20)

df_top <- df_long %>%
  filter(Species %in% top20$Species) %>%
  mutate(LogAbundance = log10(Abundance + 1))

ggplot(df_top, aes(SampleID, Species, fill=LogAbundance)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
