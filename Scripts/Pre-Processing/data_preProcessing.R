#comment these out when submitting job
rm(list=ls())
project_path=dirname(rstudioapi::getActiveDocumentContext()$path)
##################


library(glmnet)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(mltools)
library(missForestPredict)
library(caret)
library(pROC)
library(fake)
library(igraph)
library(pheatmap)
library(sharp)
library(onehot)


biomarkers <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_internal.rds"))
exposures <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_external.rds"))
ckd_case <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/outcome_definition/Outputs/output_final.rds"))

# create case-control column in exposure/biobank data, we will only use incident cases 

ckd_case$case<- as.factor(ckd_case$case)
ckd_case$prevalent_case <- as.factor(ckd_case$prevalent_case)
ckd_case$incident_case <- as.factor(ckd_case$incident_case)

table(ckd_case$case)
table(ckd_case$incident_case)
table(ckd_case$prevalent_case)

#remove prevalent cases

ckd_case_incident <- filter(ckd_case, prevalent_case == 0)

#remove prevalent cases from biomarker and exposure data
keep_ppl <- rownames(ckd_case_incident)
exposures <- exposures[keep_ppl, ]
biomarkers <- biomarkers[keep_ppl, ]

#add case_control columns to biomarker and exposure data 

exposures$case_control <- ckd_case_incident$case
biomarkers$case_control <- ckd_case_incident$case

#enforce inclusion/exclusion criteria: 

#exclusion_case_controls <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/outcome_definition/Outputs/exclusion_criteria_final.rds")
#table(exclusion_case_controls$case == 1)
#table(exclusion_case_controls$case == 1 & ckd_case$case == 1)

#need to remove 24912 individuals, 1331 who are incident CKD cases. 

#inclusion_individuals <- rownames(exclusion_case_controls)[exclusion_case_controls$case == 0]
#exposures <- exposures[inclusion_individuals,]
#biomarkers <- biomarkers[inclusion_individuals,]

#remove exposure variables with high missingness >50% 

column_missing_exp <- colMeans(is.na(exposures)) * 100
column_missing_exp <- as.data.frame(column_missing_exp) # all are <=50% therefore no need to remove 

# remove rheumatoid factor and oestradiol biomarker variables due to high missingness 

column_missing_bio <- as.data.frame(colMeans(is.na(biomarkers)) * 100)
keep_vars <- rownames(column_missing_bio)[column_missing_bio[,1] < 50 ]
biomarkers <- biomarkers[,keep_vars]

#remove individuals with high exposure missingness = 454 ppl 
individual_missing_exp <- as.data.frame(rowMeans(is.na(exposures)) * 100)
sum(individual_missing_exp$`rowMeans(is.na(exposures)) * 100` >= 50 & exposures$case_control == 0)

keep_individuals <- rownames(individual_missing_exp)[individual_missing_exp[,1] < 50 ]
exposures <- exposures[keep_individuals,]

biomarkers <- biomarkers[keep_individuals,] #carry through these removals to the biomarkers too 

#remove people with high biomarker missingness == 32,861 people, 1696 who are cases

individual_missing_bio <- as.data.frame(rowMeans(is.na(biomarkers)) * 100)
sum(individual_missing_bio$`rowMeans(is.na(biomarkers)) * 100` >= 50 & biomarkers$case_control == 0)

keep_individuals_bio <- rownames(individual_missing_bio)[individual_missing_bio[,1] < 50 ]
biomarkers <- biomarkers[keep_individuals_bio,]
exposures <- exposures[keep_individuals_bio,] # carry thorugh changes to exposure data

table(exposures$case_control)

#test train splits for exposures 

set.seed(8960)  

# Separate features and labels
Y <- as.factor(exposures$case_control)
X <- as.data.frame(exposures[, !colnames(exposures) %in% c("case_control")])

# Create index for 80/20 split
train_indices <- sample(seq_len(nrow(X)), size = 0.8 * nrow(X))

# Split the data
X_train_exp <- X[train_indices, ]
Y_train_exp <- Y[train_indices]

X_test_exp <- X[-train_indices, ]
Y_test_exp <- Y[-train_indices]

#make sure that biomarkers follow same test train split (i.e. same individuals)

Y_biomarkers <- as.factor(biomarkers$case_control)
X_biomarkers <- as.data.frame(biomarkers[, !colnames(biomarkers) %in% c("case_control")])

X_train_biomarkers <- X_biomarkers[train_indices, ]
Y_train_biomarkers <- Y_biomarkers[train_indices]

X_test_biomarkers <- X_biomarkers[-train_indices, ]
Y_test_biomarkers <- Y_biomarkers[-train_indices]

#save these for imputation using MICE

save_path <- "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/"
saveRDS(X_train_exp, file = file.path(save_path, "X_train_exp.rds"))
saveRDS(Y_train_exp, file = file.path(save_path, "Y_train_exp.rds"))

saveRDS(X_test_exp, file = file.path(save_path, "X_test_exp.rds"))
saveRDS(Y_test_exp, file = file.path(save_path, "Y_test_exp.rds"))

saveRDS(X_train_biomarkers, file = file.path(save_path, "X_train_biomarkers.rds"))
saveRDS(Y_train_biomarkers, file = file.path(save_path, "Y_train_biomarkers.rds"))

saveRDS(X_test_biomarkers, file = file.path(save_path, "X_test_biomarkers.rds"))
saveRDS(Y_test_biomarkers, file = file.path(save_path, "Y_test_biomarkers.rds"))

#imputation for exposures, impute test data based of train 
### THIS TAKES A VERY LONG TIME, PROBS WONT USE 

#imputed_exposures_train <- missForest(X_train_exp, save_models = TRUE)
#imputed_biomarkers_train <- missForest(X_train_biomarkers, save_models = TRUE)

#saveRDS(imputed_exposures_train$ximp, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/imputed_exposures_train.rds")
#saveRDS(imputed_biomarkers_train$ximp, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/imputed_biomarkers_train.rds")

#impute test data

#imputed_exposures_test <- missForestPredict(imputed_exposures_train, newdata = X_test_exp)
#imputed_biomarkers_test <- missForestPredict(imputed_biomarkers_train, newdata = X_test_biomarkers)

#saveRDS(imputed_exposures_test, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/imputed_exposures_test.rds")
#saveRDS(imputed_biomarkers_test, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/imputed_biomarkers_test.rds")

#------------------------#
#LOAD IMPUTED DATA
#------------------------#

#Biomarkers 

x_bio_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_test_biomarkers_imputed_mice.rds")
x_bio_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_biomarkers_imputed_mice.rds")

#standardise biomarker data

train_numeric_vars_bio <- x_bio_train %>% select_if(is.numeric)

Train_means_bio <- data.frame(as.list(train_numeric_vars_bio %>% apply(2, mean)))
Train_stddevs_bio <- data.frame(as.list(train_numeric_vars_bio %>% apply(2, sd)))

col_names_bio <- names(train_numeric_vars_bio)
names(Train_means_bio) <- col_names_bio
names(Train_stddevs_bio) <- col_names_bio

for (i in 1:length(col_names_bio)) {
  col <- col_names_bio[i]
  x_bio_train[, col] <- (x_bio_train[[col]] - Train_means_bio[[col]]) / Train_stddevs_bio[[col]]
  x_bio_test[, col]  <- (x_bio_test[[col]] - Train_means_bio[[col]]) / Train_stddevs_bio[[col]]
}

#save standardised data

saveRDS(x_bio_test, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_biomarker.rds")
saveRDS(x_bio_train, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds" )

#now load, scale and one hot encode imputed exposure data. also variable recoding? 

x_test_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_test_exp_imputed_mice.rds")
x_train_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_exp_imputed_mice.rds")

train_numeric_vars_exp <- x_train_exp %>% select_if(is.numeric)

Train_means_exp <- data.frame(as.list(train_numeric_vars_exp %>% apply(2, mean)))
Train_stddevs_exp <- data.frame(as.list(train_numeric_vars_exp %>% apply(2, sd)))

col_names_exp <- names(train_numeric_vars_exp)
names(Train_means_exp) <- col_names_exp
names(Train_stddevs_exp) <- col_names_exp

for (i in 1:length(col_names_exp)) {
  col <- col_names_exp[i]
  x_train_exp[, col] <- (x_train_exp[[col]] - Train_means_exp[[col]]) / Train_stddevs_exp[[col]]
  x_test_exp[, col]  <- (x_test_exp[[col]] - Train_means_exp[[col]]) / Train_stddevs_exp[[col]]
}

#one hot encoding 
#set rowname as ID colum as one hot coding messes up the rownames (i.e. the eid)

x_train_exp$eid <- as.numeric(rownames(x_train_exp))
x_test_exp$eid <- as.numeric(rownames(x_test_exp))

select_cols <- names(x_train_exp)[sapply(x_train_exp, is.factor)]
select_cols <- names(x_test_exp)[sapply(x_test_exp, is.factor)]

x_train_exp_onehot <- dummy_cols(x_train_exp, 
                        select_columns = select_cols, 
                        remove_selected_columns = TRUE, 
                        remove_first_dummy = TRUE)
x_test_exp_onehot <- dummy_cols(x_test_exp, 
                                 select_columns = select_cols, 
                                 remove_selected_columns = TRUE, 
                                 remove_first_dummy = TRUE)


##dont one hot encode this way, it keeps the reference catagory in 

#encoder <- onehot(x_train_exp, stringsAsFactors = TRUE, max_levels = Inf)
#x_train_exp_onehot <- as.data.frame(predict(encoder, x_train_exp))
#x_test_exp_onehot  <- as.data.frame(predict(encoder, x_test_exp))

saveRDS(x_test_exp_onehot, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_exposure.rds")
saveRDS(x_train_exp_onehot, "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_exposure.rds")



















