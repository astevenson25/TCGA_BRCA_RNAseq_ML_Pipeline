library(caret)
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)

# ---- load processed data ----
expr <- readRDS("data/processed/expr_filtered.rds")  
se   <- readRDS("data/processed/brca_se.rds")

# ---- extract labels, make sure they're considered factors ----
labels <- colData(se)$sample_type
labels <- factor(make.names(labels))

# ---- train/test split ----
set.seed(123)
train_index <- createDataPartition(labels, p=0.8, list = FALSE)

train_expr_initial <- expr[, train_index]
test_expr_initial <- expr[, -train_index]

train_labels <- labels[train_index]
test_labels  <- labels[-train_index]

# ---- feature selection (training set only) ----

# variance filtering, top 500 genes
gene_variances <- apply(train_expr_initial, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(gene_variances, decreasing = TRUE))[1:500]

train_expr_var <- train_expr_initial[top_var_genes,]

train_normal_idx <- which(train_labels == "Solid.Tissue.Normal")
train_tumor_idx <- which(train_labels == "Primary.Tumor")

# Ensure column names are characters
colnames(train_expr_var) <- as.character(colnames(train_expr_var))

# differential expression on training set
DEGs_train <- TCGAanalyze_DEA(
  mat1 = train_expr_var[, train_normal_idx],
  mat2 = train_expr_var[, train_tumor_idx],
  pipeline = "limma",
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT",
  ClinicalDF = data.frame()
)

selected_genes <- rownames(DEGs_train)

# ---- subset train and test using selected genes ----
train_expr <- train_expr_initial[selected_genes,]
test_expr <- test_expr_initial[selected_genes,]

# ----  CV folds (fixed across models) ----
set.seed(123)
folds <- createFolds(train_labels, k = 5, returnTrain = TRUE)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  index = folds,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

dir.create("data/model_inputs", recursive = TRUE, showWarnings = FALSE)

saveRDS(train_expr,  "data/model_inputs/train_expr.rds")
saveRDS(test_expr,   "data/model_inputs/test_expr.rds")
saveRDS(train_labels,"data/model_inputs/train_labels.rds")
saveRDS(test_labels, "data/model_inputs/test_labels.rds")
saveRDS(DEGs_train, "data/model_inputs/DEGs_train.rds")
saveRDS(selected_genes, "data/model_inputs/selected_genes.rds")
saveRDS(ctrl, "data/model_inputs/ctrl.rds")



