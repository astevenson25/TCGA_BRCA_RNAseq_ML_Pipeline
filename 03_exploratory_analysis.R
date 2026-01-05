library(tidyverse)
library(pheatmap)

# ---- load in saved data from preprocessing script ----
expr <- readRDS("data/processed/expr_filtered.rds")
se   <- readRDS("data/processed/brca_se.rds")

# ---- extract clinical metadata ----
clin <- as.data.frame(colData(se))

# ---- keep only samples present in expression matrix ----
clin <-clin[match(colnames(expr), clin$barcode),]

# ---- select key clinical variables ----
clin <- clin %>%
  mutate(
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    days_to_death = as.numeric(days_to_death),
    sample_type = shortLetterCode
  )

# ---- select highly variable genes for EDA ----
gene_var <- apply(expr, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:500]

expr_var <- expr[top_genes,]

# ---- PCA analysis ----
expr_t <- t(expr_var)
pca <- prcomp(expr_t, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample_type = clin$sample_type,
  days_to_death = clin$days_to_death
)

# ---- save PCA plot ----
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
png("results/plots/pca_top500_variable_genes.png", width = 900, height = 600)
plot(
  pca_df$PC1, pca_df$PC2,
  col = ifelse(pca_df$sample_type == "TP", "red", "blue"),
  pch = 16,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA of Top 500 Most Variable Genes"
)
legend("topright", legend = c("Tumor", "Normal"),
       col = c("red", "blue"), pch = 16)
dev.off()


# ---- cumulative variance explained ----
var_explained <- cumsum(pca$sdev^2 / sum(pca$sdev^2))

png("results/plots/pca_cumulative_variance.png", width = 800, height = 600)
plot(var_explained, type = "b",
     xlab = "Principal Component",
     ylab = "Cumulative Variance Explained",
     main = "PCA cumulative variance")
abline(h = 0.75, col = "red", lty = 2)
dev.off()

# ---- gene-survival correlation (Spearman) ----
surv <- clin$days_to_death
cor_res <- apply(expr_var, 1, function(g){
  if (sd(g, na.rm = TRUE) == 0) return(c(rho = NA, pval = NA))
  test <- cor.test(g, surv, method = "spearman", use = "pairwise.complete.obs")
  c(rho = test$estimate, pval = test$p.value)
})

cor_df <- as.data.frame(t(cor_res))
cor_df$padj <- p.adjust(cor_df$pval, method = "fdr")

# ---- create and save heatmap of top correlated genes ----
top_cor_genes <- rownames(cor_df)[order(cor_df$padj)][1:50]

rho_mat <- matrix(cor_df[top_cor_genes, "rho.rho"], ncol = 1)
rownames(rho_mat) <- top_cor_genes
colnames(rho_mat) <- "days_to_death"

png("results/plots/gene_survival_correlation_heatmap.png", width = 600, height = 900)
pheatmap(
  rho_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression vs Survival (Spearman rho)"
)
dev.off()


# ---- save correlation results ----
saveRDS(cor_df, "results/deg_tables/gene_survival_correlations.rds")

