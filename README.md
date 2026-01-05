# TCGA-BRCA RNA-seq Analysis Pipeline
This repository contains a pipeline for processing TCGA Breast Cancer (BRCA) RNA-seq data, performing exploratory data analysis, feature selection, and training/testing multiple machine learning models.


# Table of Contents
1. Overview
2. Pipeline Steps
3. Software Requirements
4. Directory Structure
5. Usage
6. Outputs
7. Notes

# 1. Overview

This pipeline processes TCGA-BRCA RNA-seq data to generate machine learning models that distinguish primary breast cancer tumor tissue from normal breast tissue.

General steps include:
- downloading data from TCGA
- preprocessing and filtering of protein-coding genes
- exploratory data analysis
- feature selection
- model training and evaluation (logistic regression, random forest, radial SVM, and neural network)
- visualization of model performance results

# 2. Pipeline Steps

01_download_data.R
- Purpose: Downloads TCGA BRCA RNA-seq data as processed unstranded FPKM values (as provided by TCGA for consistency and computational efficiency) for 100 primary tumor and 100 normal tissue samples (downsampled for memory constraints), along with clinical metadata. Stores results in a SummarizedExperiment object.
- Input: 
    - TCGA project ID (TCGA-BRCA)
- Output: 
    - data/processed/brca_se.rds
    - GDCdata directory

02_preprocessing.R
- Purpose: Filters protein-coding genes, removes genes with low correlation across samples (cor.cut = 0.6), and removes the lowest 25% of genes by expression using FPKM values. Differential expression analysis is performed using the 'TCGAanalyze_DEA' function. Boxplots of raw and filtered expression data are generated to assess normalization and filtering outcomes.
- Input:
    - data/processed/brca_se.rds
- Output:
    - data/processed/expr_filtered.rds
    - results/deg_tables/degs_full.rds
    - results/plots/filtered_expression_boxplot.png
    - results/plots/raw_expression_boxplot.png

03_exploratory_analysis.R
- Purpose: Perform exploratory PCA and gene-survival correlations. Outputs PCA of top 500 genes with the most variance, a cumulative variance plot, a gene-survival correlation heatmap, and saves the table of gene-survival correlations.
- Input:
    - data/processed/expr_filtered.rds
    - data/processed/brca_se.rds
- Output:
    - results/plots/pca_top500_variable_genes.png
    - results/plots/pca_cumulative_variance.png
    - results/plots/gene_survival_correlation_heatmap.png
    - results/deg_tables/gene_survival_correlations.rds

04_feature_selection.R
- Purpose: Splits data into training and testing sets and performs feature selection using only the training data to avoid data leakage. Feature selection is based on variance filtering followed by differential expression analysis. The selected genes are then used to subset both training and testing datasets.

- Input:
    - data/processed/expr_filtered.rds
    - data/processed/brca_se.rds
- Output:
    - data/model_inputs/train_expr.rds
    - data/model_inputs/test_expr.rds
    - data/model_inputs/train_labels.rds
    - data/model_inputs/test_labels.rds
    - data/model_inputs/selected_genes.rds
    - data/model_inputs/DEGs_train.rds
    - data/model_inputs/ctrl.rds

05_model_training.R
- Purpose: Trains classification models (logistic regression, random forest, radial SVM, and neural network) using 5-fold cross validation and tests model performance. A random seed is set prior to training to ensure reproducibility of results. Performance metrics are output for each model.
- Input:
    - data/model_inputs/train_expr.rds
    - data/model_inputs/test_expr.rds
    - data/model_inputs/train_labels.rds
    - data/model_inputs/test_labels.rds
    - data/model_inputs/ctrl.rds
- Output:
    - results/models/log_model.rds
    - results/models/rf_model.rds
    - results/models/svm_model.rds
    - results/models/nn_model.rds
    - results/models/test_predicted_probabilities.rds
    - results/models/auc_summary.rds
    - results/models/confusion_matrices.rds

06_visualization.R
-  Purpose: Visualize and compare model performance. ROC curves, confusion matrices, and a bar graph summarizing confusion matrix metrics across models are produced.
- Input:
    - results/models/test_predicted_probabilities.rds
    - results/models/auc_summary.rds
    - results/models/confusion_matrices.rds
- Output:
    - results/models/sensitivity_specificity.rds
    - results/plots/roc_curves.png
    - results/plots/confusion_matrices.png
    - results/plots/confusion_matrix_comparison.png
    - results/plots/sensitivity_specificity.png

# Software Requirements
- R >= 4.2
- Packages:
    - TCGAbiolinks
    - tidyverse
    - EDASeq
    - SummarizedExperiment
    - caret
    - pROC
    - pheatmap
    - glmnet
    - randomForest
    - nnet
    - ggplot2
    - gridExtra
    - here

You can install packages using:
```r
BiocManager::install(c("TCGAbiolinks", "EDASeq", "SummarizedExperiment"))

install.packages(c("tidyverse", "caret", "pROC", "pheatmap", "glmnet", "randomForest", "nnet", "ggplot2", "gridExtra", "here", "BiocManager"))

```
# Directory Structure

```
project_root/
├── data/
│   ├── processed/      # Processed RDS objects
│   └── model_inputs/   # Train/test sets for ML
├── results/
│   ├── deg_tables/     # Differential expression results
│   ├── plots/          # PCA, boxplots, heatmaps, ROC, CM plots
│   └── models/         # Trained ML models
├── GDCdata/
│   └── TCGA-BRCA/      # TCGA-BRCA data
├── 01_download_data.R
├── 02_preprocessing.R
├── 03_exploratory_analysis.R
├── 04_feature_selection.R
├── 05_model_training.R
└── 06_visualization.R
```

# Usage

In R, run the pipeline scripts in the following order:
```r

source("01_download_data.R")
source("02_preprocessing.R")
source("03_exploratory_analysis.R")
source("04_feature_selection.R")
source("05_model_training.R")
source("06_visualization.R")

```
# Outputs
 - Processed data: brca_se.rds, expr_filtered.rds
 - DEGs tables: degs_full.rds 
 - EDA plots: PCA, boxplots, gene-survival heatmaps
 - ML inputs: Train/test matrices, selected genes, CV folds
 - Trained models: log_model.rds, rf_model.rds, svm_model.rds, nn_model.rds
 - Model evaluation: ROC curves, confusion matrices, AUC summary

 # Notes
 - Downsampling is applied to the original dataset to include 100 tumor and 100 normal tissue samples for memory constraints and sample balance. Consider removing or increasing this for full dataset analysis.
 - Feature selection is currently based on the top 500 most variable differentially expressed genes from the training set.
 - For reproducible results, a random seed is set during model training and testing. This ensures consistent train/test splits and model outputs.
  - If any changes are made to sample size, ensure factor levels for sample labels are consistent across train/test sets.
 - Make sure working directory is set appropriately.

