library(EDASeq)
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

# ---- load in saved data from download_data script ----
data <- readRDS("data/processed/brca_se.rds")

# ---- create folders if they don’t exist ----
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/deg_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# ---- subset to protein coding genes only ----
SECoding <- data[rowData(data)$gene_type == "protein_coding", ]

### ---- filtering ---- ###

# ---- data preprocessing for protein coding genes ----
# uses FPKM values from unstranded sequencing workflow
# filters out genes with low correlation across samples (cutoff >= 0.6)
dataPrep_Coding <- TCGAanalyze_Preprocessing(
  object = SECoding, 
  cor.cut = 0.6,  
  datatype = "fpkm_unstrand"
)

# method quantile and cutoff = 0.25 removes lowest 25% of genes by expression
dataFilt_Coding <- TCGAanalyze_Filtering(
  tabDF = dataPrep_Coding, 
  method = "quantile", 
  qnt.cut =  0.25
)   

# ---- identify sample groups ----
sample_types <- colData(data)$shortLetterCode

dataNormal_Target  <- colnames(dataFilt_Coding)[sample_types == "NT"]
dataPrimary_Target <- colnames(dataFilt_Coding)[sample_types == "TP"]


# ---- differential expression analysis ----
DEGsCoding <- TCGAanalyze_DEA(mat1 = dataFilt_Coding[,dataNormal_Target],
                              mat2 = dataFilt_Coding[,dataPrimary_Target],
                              pipeline="limma",
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT", 
                              ClinicalDF = data.frame()
)

# ---- save filtered expression matrix ----
saveRDS(dataFilt_Coding, "data/processed/expr_filtered.rds")

# ---- save DEGs table ----
saveRDS(DEGsCoding, "results/deg_tables/degs_full.rds")

# ---- save boxplots of raw and filtered expression data ----
png("results/plots/raw_expression_boxplot.png", width = 900, height = 600)
boxplot(assay(SECoding),
        outline = FALSE,
        las = 2,
        main = "Raw FPKM expression (protein-coding genes)")
dev.off()

png("results/plots/filtered_expression_boxplot.png", width = 900, height = 600)
boxplot(dataFilt_Coding,
        outline = FALSE,
        las = 2,
        main = "Filtered FPKM expression (protein-coding genes)")
dev.off()



