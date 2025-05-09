
library(mice)
library(parallel)
library(doParallel)

#load exposure data

# Define path
load_path <- "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/"

# Load exposure train/test data
X_train_exp <- readRDS(file.path(load_path, "X_train_noimpute_exp.rds"))
X_test_exp <- readRDS(file.path(load_path, "X_test_noimpute_exp.rds"))

#--------------------------------------------------
# Set up parallel processing
#--------------------------------------------------
n_cores <- 29  # Leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)

#--------------------------------------------------
# Imputation for Exposure data
#--------------------------------------------------
# Create predictor matrix for training exposure data
pred_matrix_exp <- make.predictorMatrix(X_train_exp)

# Run MICE on training exposure data
mice_train_exp <- mice(X_train_exp, 
                       m = 5,               # Number of multiple imputations
                       maxit = 5,           # Number of iterations
                       predictorMatrix = pred_matrix_exp,
                       printFlag = TRUE,    # Set to TRUE to monitor progress
                       n.core = n_cores,    # Use multiple cores
                       n.imp.core = 1)      # Imputations per core

# Extract complete imputed training dataset
X_train_exp_imputed <- complete(mice_train_exp, 1)

saveRDS(X_train_exp_imputed, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_exp_imputed_mice.rds")

# Apply same imputation model to test exposure data
mice_test_exp <- mice.mids(mice_train_exp, newdata = X_test_exp)
X_test_exp_imputed <- complete(mice_test_exp, 1)

saveRDS(X_test_exp_imputed, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_test_exp_imputed_mice.rds")
