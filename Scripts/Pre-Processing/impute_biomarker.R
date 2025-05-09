library(mice)
library(parallel)
library(doParallel)

#--------------------------------------------------
# Load biomarker data
#--------------------------------------------------
# Define path
load_path <- "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/"

# Load biomarker train/test data
X_train_biomarkers <- readRDS(file.path(load_path, "X_train_noimpute_biomarkers.rds"))
X_test_biomarkers <- readRDS(file.path(load_path, "X_test_noimpute_biomarkers.rds"))

#--------------------------------------------------
# Set up parallel processing
#--------------------------------------------------
n_cores <- 29
cl <- makeCluster(n_cores)
registerDoParallel(cl)

#--------------------------------------------------
# Imputation for Biomarker data
#--------------------------------------------------
# Create predictor matrix for training biomarker data
pred_matrix_biomarkers <- make.predictorMatrix(X_train_biomarkers)

# Run MICE on training biomarker data
mice_train_biomarkers <- mice(X_train_biomarkers, 
                              m = 5,               # Number of multiple imputations
                              maxit = 5,           # Number of iterations
                              predictorMatrix = pred_matrix_biomarkers,
                              printFlag = TRUE,    # Set to TRUE to monitor progress
                              n.core = n_cores,    # Use multiple cores
                              n.imp.core = 1)      # Imputations per core

# Extract complete imputed training dataset
X_train_biomarkers_imputed <- complete(mice_train_biomarkers, 1)

saveRDS(X_train_biomarkers_imputed, file.path(load_path, "X_train_biomarkers_imputed_mice.rds"))

# Apply same imputation model to test biomarker data
mice_test_biomarkers <- mice.mids(mice_train_biomarkers, newdata = X_test_biomarkers)
X_test_biomarkers_imputed <- complete(mice_test_biomarkers, 1)

saveRDS(X_test_biomarkers_imputed, file.path(load_path, "X_test_biomarkers_imputed_mice.rds"))
