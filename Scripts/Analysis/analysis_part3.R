library(glmnet)
library(caret)
library(pROC)

###------------------------- Step 9---------------------------------

# Load all selected variables
selected_vars1 <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars1.csv")
selected_vars_HDL <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_HDL.csv")
selected_vars_RBC <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_RBC.csv")

# Load datasets
Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")
Y_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_test.rds")

table(Y_train)
table(Y_test)

# Combine all selected exposures
all_selected_exp <- c(as.character(selected_vars1$variable),
                      as.character(selected_vars_HDL$variable),
                      as.character(selected_vars_RBC$variable))
all_selected_exp <- unique(all_selected_exp)
list(all_selected_exp)

# Add confounders
confounders <- c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black")
selected_biomarkers <- c("HDL_cholesterol", "nucleated_red_blood_cell_count")

# Load original dataset
X_train_exposure <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_exposure.rds")
X_train_biomarker <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")

# train/validation split for the origin datasets
set.seed(2025)

# split
n <- nrow(X_train_exposure)
train_prop <- 0.625
val_prop <- 0.375
train_size <- round(train_prop * n)

all_indices <- sample(1:n)
train_indices <- all_indices[1:train_size]
val_indices <- all_indices[(train_size + 1):n]

# train sets
X_train50_exp <- X_train_exposure[train_indices, ]
X_train50_bio <- X_train_biomarker[train_indices, ]
Y_train50 <- Y_train[train_indices]
table(Y_train50)

# validation sets
X_val30_exp <- X_train_exposure[val_indices, ]
X_val30_bio <- X_train_biomarker[val_indices, ]
Y_val30 <- Y_train[val_indices]
table(Y_val30)

# save
saveRDS(X_train50_exp, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train50_exp.rds")
saveRDS(X_train50_bio, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train50_bio.rds")
saveRDS(Y_train50, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train50.rds")

saveRDS(X_val30_exp, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_val30_exp.rds")
saveRDS(X_val30_bio, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_val30_bio.rds")
saveRDS(Y_val30, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_val30.rds")

# Prepare train and test sets
X_test_exposure <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_exposure.rds")
X_test_biomarker <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_biomarker.rds")
Y_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_test.rds")

X_train_selected_exp <- X_train50_exp[, all_selected_exp]
X_test_selected_exp <- X_test_exposure[, all_selected_exp]

X_train_conf <- X_train50_exp[, confounders]
X_test_conf <- X_test_exposure[, confounders]

X_train_selected_bio <- X_train50_bio[, selected_biomarkers]
X_test_selected_bio <- X_test_biomarker[, selected_biomarkers]

# Combine
X_train_complete <- cbind(X_train_selected_exp, X_train_conf, X_train_selected_bio)
X_test_complete <- cbind(X_test_selected_exp, X_test_conf, X_test_selected_bio)

saveRDS(X_train_complete, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_complete.rds")
saveRDS(X_test_complete, file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_test_complete.rds")

# transform 
X_train_complete <- as.matrix(X_train_complete)
X_test_complete <- as.matrix(X_test_complete)
Y_train50 <- as.numeric(Y_train50)
Y_test <- as.numeric(Y_test)
Y_train50 <- ifelse(Y_train50 == 2, 1, 0)
Y_test <- ifelse(Y_test == 2, 1, 0)
table(Y_train50)
table(Y_test)

##------ Lasso--------
set.seed(202)
lasso_complete_cv <- cv.glmnet( x = X_train_complete, y = Y_train50, alpha = 1, family = "binomial")
lasso_complete <- glmnet(X_train_complete, Y_train50, alpha = 1, lambda = lasso_complete_cv$lambda.1se, family = "binomial")

# extract nonzero vars
table(coef(lasso_complete, s = "lasso_complete$lambda.1se")[-1] != 0)

# extract nonzero coefficients
lasso_coefs_complete <- as.matrix(coef(lasso_complete, s = "lambda.1se"))
lasso_coefs_nonzero <- lasso_coefs_complete[lasso_coefs_complete[,1] != 0, , drop = FALSE]  
lasso_coefs_nonzero <- lasso_coefs_nonzero[rownames(lasso_coefs_nonzero) != "(Intercept)", , drop = FALSE]  

# sort coefficients by absolute value
lasso_complete_coef_sorted <- lasso_coefs_nonzero[order(abs(lasso_coefs_nonzero[,1]), decreasing = TRUE), , drop = FALSE]
lasso_complete_coef_sorted <- data.frame( Variable = rownames(lasso_complete_coef_sorted), Coefficient = as.numeric(lasso_complete_coef_sorted))

# show coefficients
print(lasso_complete_coef_sorted)

# save the coefficient results
write.csv(as.data.frame(lasso_complete_coef_sorted), 
          file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/lasso_coefficients_sorted.csv")

##--------Elastic-net--------

# CV error for a given alpha
cv_enet <- function(alpha) {
  cv_model <- cv.glmnet(x = X_train_complete, y = Y_train50, alpha = alpha, family = "binomial", nfolds = 10)
  cv_error <- cv_model$cvm[which.min(abs(cv_model$lambda - cv_model$lambda.1se))]
  return(cv_error)
}

# optimise alpha
set.seed(2025)
alpha_opt <- optimise(cv_enet, c(0, 1))
optimal_alpha <- alpha_opt$minimum
cat("Optimal alpha:", optimal_alpha, "\n")

# elastic net
enet_cv <- cv.glmnet(x = X_train_complete, y = Y_train50, alpha = optimal_alpha, family = "binomial", nfolds = 10)
enet <- glmnet(X_train_complete, Y_train50, alpha = optimal_alpha, lambda = enet_cv$lambda.1se, family = "binomial")

# extract coefficients
enet_coefs <- coef(enet_cv, s = "lambda.1se")
enet_betas <- enet_coefs[-1]  
names(enet_betas) <- rownames(enet_coefs)[-1]

# non-zero coefficients
enet_coefs_nonzero <- enet_betas[enet_betas != 0]
enet_coefs_table <- data.frame(Variable = names(enet_coefs_nonzero), 
                               Coefficient = as.numeric(enet_coefs_nonzero))
enet_coefs_sorted <- enet_coefs_table[order(abs(enet_coefs_table$Coefficient), decreasing = TRUE), ]

print(enet_coefs_sorted)

# save the coefficient results
write.csv(as.data.frame(enet_coefs_sorted), 
          file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/elastic_coefficients_sorted.csv")

#-----------AUC--------------
# prediction
lasso_pred_prob <- predict(lasso_complete, newx = X_test_complete, s = "lambda.1se", type = "response")

# predicted probabilities to class
lasso_pred_class <- as.factor(ifelse(lasso_pred_prob > 0.5, 1, 0))

# Lasso AUC
roc_lasso <- roc(Y_test, as.numeric(lasso_pred_prob))
auc_lasso <- auc(roc_lasso)
cat("Lasso AUC:", auc_lasso, "\n")

# Confusion matrix
Y_test <- factor(Y_test, levels = c(0, 1))  
confusionMatrix(lasso_pred_class, Y_test, positive = "1")

# prediction
enet_pred_prob <- predict(enet, newx = X_test_complete, s = "lambda.1se", type = "response")
enet_pred_class <- as.factor(ifelse(enet_pred_prob > 0.5, 1, 0))

# Elastic-net AUC
roc_enet <- roc(Y_test, as.numeric(enet_pred_prob))
auc_enet <- auc(roc_enet)
cat("Elastic-net AUC:", auc_enet, "\n")

# Confusion matrix
confusionMatrix(enet_pred_class, Y_test, positive = "1")


