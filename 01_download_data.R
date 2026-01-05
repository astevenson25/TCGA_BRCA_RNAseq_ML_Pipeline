library(TCGAbiolinks)
library(tidyverse)
library(here)

# ---- obtain project summary information----
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

# ---- initial query to get metadata for the case ----
query_meta <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

CaseInfo <- getResults(query_meta)
as_tibble(head(CaseInfo))

# ---- select tumor and normal samples ----

# primary tumor tissue
dataPrimary_Target <- TCGAquery_SampleTypes(
  barcode = CaseInfo$cases, 
  typesample = "TP"
)

# normal tissue
dataNormal_Target <- TCGAquery_SampleTypes(
  barcode = CaseInfo$cases,
  typesample = "NT"
) 

# ---- downsample due to memory constraints ----
dataPrimary_Target <- dataPrimary_Target[1:100]
dataNormal_Target <- dataNormal_Target[1:100]

# ---- query expression data for the selected samples ----
TargetSamples <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(dataPrimary_Target, dataNormal_Target)
)


# ---- download and prepare data ----
message("Downloading TGCA BRCA RNAseq data. . .")
GDCdownload(TargetSamples) 

# ---- create summarized experiment GDS object and directory to store results ----
message("Preparing summarized experiment object. . .")
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
se <- GDCprepare(TargetSamples)

# ---- save GDS object output ----
message("Saving processed object. . .")
saveRDS(se, here("data", "processed", "brca_se.rds"))


