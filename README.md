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
