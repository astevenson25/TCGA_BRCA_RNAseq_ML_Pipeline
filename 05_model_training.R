library(caret)
library(dplyr)
library(pROC)

set.seed(123)

# ---- load model inputs from feature selection ----
train_expr   <- readRDS("data/model_inputs/train_expr.rds")
test_expr    <- readRDS("data/model_inputs/test_expr.rds")
train_labels <- readRDS("data/model_inputs/train_labels.rds")
test_labels  <- readRDS("data/model_inputs/test_labels.rds")
ctrl <- readRDS("data/model_inputs/ctrl.rds")

# ---- checks that factor levels in train set match test set ----
stopifnot(identical(levels(train_labels), levels(test_labels)))

# ---- ensure tumor tissue is positive class for ROC/AUC ----
train_labels <- factor(train_labels, levels = c("Solid.Tissue.Normal", "Primary.Tumor"))
test_labels  <- factor(test_labels,  levels = c("Solid.Tissue.Normal", "Primary.Tumor"))
positive_class <- levels(train_labels)[2]  # Primary.Tumor

# ---- transpose expression matrices (samples x genes) ----
X_train <- t(train_expr)
X_test  <- t(test_expr)

# ---- train models ----

# logistic regression (elastic net)
set.seed(123)
log_model <- train(
  x = X_train,
  y = train_labels,
  method = "glmnet",
  trControl = ctrl,
  metric = "ROC"
)

# random forest
set.seed(123)
rf_model <- train(
  x = X_train,
  y = train_labels,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC"
)

# support vector machine (radial)
set.seed(123)
svm_model <- train(
  x = X_train,
  y = train_labels,
  method = "svmRadial",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC"
)

# neural network
set.seed(123)
nn_model <- train(
  x = X_train,
  y = train_labels,
  method = "nnet",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC",
  trace = FALSE
)

### ---- test set evaluation ---- ####

# ---- class predictions ----
pred_log <- predict(log_model, X_test)
pred_rf <- predict(rf_model, X_test)
pred_svm <- predict(svm_model, X_test)
pred_nn <- predict(nn_model, X_test)

# ---- confusion matrices ----
cm_log <- confusionMatrix(pred_log, test_labels)
cm_rf <- confusionMatrix(pred_rf, test_labels)
cm_svm <- confusionMatrix(pred_svm, test_labels)
cm_nn <- confusionMatrix(pred_nn, test_labels)

# ---- predicted probabilities ----
pred_probs <- data.frame(
  logistic = predict(log_model, X_test, type = "prob")[, positive_class],
  rf = predict(rf_model, X_test, type = "prob")[, positive_class],
  svm = predict(svm_model, X_test, type = "prob")[, positive_class],
  nn = predict(nn_model, X_test, type = "prob")[, positive_class],
  label = test_labels
)

# ---- ROC/AUC ----
roc_log <- roc(test_labels, pred_probs$logistic)
roc_rf <- roc(test_labels, pred_probs$rf)
roc_svm <- roc(test_labels, pred_probs$svm)
roc_nn <- roc(test_labels, pred_probs$nn)

auc_table <- data.frame(
  Model = c("Logistic", "RandomForest", "SVM", "NeuralNetwork"),
  AUC = c(
    auc(roc_log),
    auc(roc_rf),
    auc(roc_svm),
    auc(roc_nn)
  )
)

### ---- save outputs for visualization script ---- ####
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(log_model, "results/models/log_model.rds")
saveRDS(rf_model,  "results/models/rf_model.rds")
saveRDS(svm_model, "results/models/svm_model.rds")
saveRDS(nn_model,  "results/models/nn_model.rds")

saveRDS(pred_probs, "results/models/test_predicted_probabilities.rds")
saveRDS(auc_table,  "results/models/auc_summary.rds")

saveRDS(
  list(
    Logistic = cm_log,
    RandomForest = cm_rf,
    SVM = cm_svm,
    NeuralNetwork = cm_nn
  ),
  "results/models/confusion_matrices.rds"
)