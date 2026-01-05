library(ggplot2)
library(dplyr)
library(pROC)


# ---- load in model outputs ----
pred_probs <- readRDS("results/models/test_predicted_probabilities.rds")
auc_table  <- readRDS("results/models/auc_summary.rds")
cms        <- readRDS("results/models/confusion_matrices.rds")

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# ---- ensure label factor order ----
pred_probs$label <- factor(pred_probs$label, levels = c("Solid.Tissue.Normal", "Primary.Tumor"))

# ---- create model color palette ----
model_colors <- c(
  Logistic       = "steelblue",
  RandomForest  = "darkgreen",
  SVM            = "firebrick",
  NeuralNetwork  = "purple"
)


# ---- create and save ROC curves ----

roc_log <- roc(pred_probs$label, pred_probs$logistic)
roc_rf <- roc(pred_probs$label, pred_probs$rf)
roc_svm <- roc(pred_probs$label, pred_probs$svm)
roc_nn <- roc(pred_probs$label, pred_probs$nn)

png("results/plots/roc_curves.png", width = 800, height = 800)

plot(
  roc_log,
  col = model_colors["Logistic"],
  lwd = 2,
  main = "ROC Curves for Tumor vs Normal Classification",
  asp = 1
)

plot(roc_rf,  col = model_colors["RandomForest"], lwd = 2, add = TRUE)
plot(roc_svm, col = model_colors["SVM"],          lwd = 2, add = TRUE)
plot(roc_nn,  col = model_colors["NeuralNetwork"],lwd = 2, add = TRUE)


legend(
  "bottomright",
  legend = paste0(
    auc_table$Model, " (AUC = ", round(auc_table$AUC, 3), ")"
  ),
  col = c("steelblue", "darkgreen", "firebrick", "purple"),
  lwd = 2,
  bty = "n"
)

dev.off()

# ---- create and save onfusion matrix heatmaps ----

# define plot_cm function
plot_cm <- function(cm, model_name) {
  df <- as.data.frame(cm$table)
  ggplot(df, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "white", size = 5) +
    scale_fill_gradient(
      low = "grey40",
      high = model_colors[model_name]
    ) +
    labs(title = model_name) +
    theme_minimal(base_size = 14)
}


png("results/plots/confusion_matrices.png", width = 1200, height = 1000)

gridExtra::grid.arrange(
  plot_cm(cms$Logistic, "Logistic"),
  plot_cm(cms$RandomForest, "RandomForest"),
  plot_cm(cms$SVM, "SVM"),
  plot_cm(cms$NeuralNetwork, "NeuralNetwork"),
  ncol = 2
)
dev.off()

# ---- create and save confusion matrix bar chart ----

extract_cm_counts <- function(cm, model_name) {
  tibble(
    Model = model_name,
    Metric = c("TP", "FP", "TN", "FN"),
    Value = c(
      cm$table["Primary.Tumor", "Primary.Tumor"],
      cm$table["Primary.Tumor", "Solid.Tissue.Normal"],
      cm$table["Solid.Tissue.Normal", "Solid.Tissue.Normal"],
      cm$table["Solid.Tissue.Normal", "Primary.Tumor"]
    )
  )
}

cm_bar <- bind_rows(
  extract_cm_counts(cms$Logistic, "Logistic"),
  extract_cm_counts(cms$RandomForest, "RandomForest"),
  extract_cm_counts(cms$SVM, "SVM"),
  extract_cm_counts(cms$NeuralNetwork, "NeuralNetwork")
)

png("results/plots/confusion_matrix_comparison.png", width = 1000, height = 700)

print(
  ggplot(cm_bar, aes(Model, Value, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~Metric) +
    scale_fill_manual(values = model_colors) +
    theme_minimal(base_size = 20) +
    labs(
      title = "Confusion Matrix Breakdown Across Models",
      y = "Sample Count",
      x = ""
    )
)  

dev.off()

# ---- extract and visualize sensitivity/specificity ----
extract_metrics <- function(cm) {
  data.frame(
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"]
  )
}

metric_df <- bind_rows(
  Logistic       = extract_metrics(cms$Logistic),
  RandomForest   = extract_metrics(cms$RandomForest),
  SVM            = extract_metrics(cms$SVM),
  NeuralNetwork  = extract_metrics(cms$NeuralNetwork),
  .id = "Model"
)

rownames(metric_df) <- 1:nrow(metric_df)

saveRDS(metric_df, "results/models/sensitivity_specificity.rds")

# convert to long format for faceting
metric_long <- metric_df %>%
  pivot_longer(
    cols = c(Sensitivity, Specificity),
    names_to = "Metric",
    values_to = "Value"
)

png("results/plots/sensitivity_specificity.png", width = 1000, height = 700)

ggplot(metric_long, aes(x = Model, y = Value, fill = Model)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Metric) +
  scale_fill_manual(values = model_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Model Sensitivity and Specificity",
    y = "Metric Value",
    x = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  "results/plots/sensitivity_specificity.png",
  width = 8,
  height = 4,
  dpi = 300
)




