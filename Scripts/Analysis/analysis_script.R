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
#library(onehot)

#load data 

X_train_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
X_test_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_biomarker.rds")

X_train_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_exposure.rds")
X_test_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_exposure.rds")

#nb y vectors are the same for exposures and biomarkers

Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")
Y_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_test.rds")

# remove age and sex and ethnicity from exposures to prevent leakage, these will be added back as confounders later 

#X_train_exp <- X_train_exp[, !(colnames(X_train_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]
#X_test_exp  <- X_test_exp[, !(colnames(X_test_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]

#add eid to biomarker columns 

rownames(X_train_exp) <- rownames(X_train_bio)
rownames(X_test_exp) <- rownames(X_test_bio)

X_train_bio$eid <- as.character(rownames(X_train_bio))
X_train_exp$eid <- as.character(rownames(X_train_exp))

X_test_bio$eid <- as.character(rownames(X_test_bio))
X_test_exp$eid <- as.character(rownames(X_test_exp))

# Merge datasets by 'eid' column
X_train_all <- merge(X_train_bio, X_train_exp, by = "eid")
X_train_all <- X_train_all[, !colnames(X_train_all) %in% c("eid")]

X_test_all <- merge(X_test_bio, X_test_exp, by = "eid")
X_test_all <- X_test_all[, !colnames(X_test_all) %in% c("eid")]

X_train_all_full <- X_train_all

## split train_all into train_all 50% and validation 30% of the dataset 
set.seed(2025)

n <- nrow(X_train_all_full)
train_prop <- 0.625
val_prop <- 0.375
train_size <- round(train_prop * n)


# Step 3: Sample indices
all_indices <- sample(1:n)
train_indices <- all_indices[1:train_size]
val_indices <- all_indices[(train_size + 1):n]

# Step 4: Create data splits
X_train_all <- X_train_all_full[train_indices, ]
X_val_all   <- X_train_all_full[val_indices, ]
#saveRDS(as.data.frame(X_val_all), file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_val_all.rds")


###################################

# format data for analysis 

# Convert X data frames to matrices
X_train_bio <- as.matrix(X_train_bio)
X_test_bio <- as.matrix(X_test_bio)

X_train_exp <- as.matrix(X_train_exp)
X_test_exp <- as.matrix(X_test_exp)

X_train_all <- as.matrix(X_train_all)
X_test_all <- as.matrix(X_test_all)

X_val_all <- as.matrix(X_val_all)

# Convert Y to numeric vectors
Y_train <- as.numeric(as.character(Y_train))
Y_test <- as.numeric(as.character(Y_test))

##############LASSO stability selection for direct model (i.e. using both exposures and outcome joint matrix) ########

set.seed(2025)

X_conf <- X_train_all[, colnames(X_train_all) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black")]
X_train_all  <- X_train_all[, !(colnames(X_train_all) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]

model1 <- VariableSelection(
  xdata = cbind(X_train_all,X_conf),
  ydata = Y_train[train_indices],
  family = "binomial",
  n_cat = 3,
  pi_list = seq(0.5, 0.99, by = 0.01), 
  penalty.factor = c(rep(1,
                         ncol(X_train_all)), rep(0, ncol(X_conf)))
)

CalibrationPlot(model1)
#png("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/calibration_plot_model1.png",
    #width = 1000, height = 800, res = 150)
#par(mar = c(5, 5, 4, 2))  # optional: adjust margins if needed
#CalibrationPlot(model1)
#dev.off()

#find the optimal pi threshold (argmax of stability score)
optimal_params1 <- Argmax(model1)
print(optimal_params1) #optimal lambda:0.0004920507 && optimal pi: 0.78

#selection proportions of all variables (not just the optimal ones)
selprop1 <- SelectionProportions(model1)
selprop1_sorted <- sort(selprop1, decreasing = TRUE)
print(selprop1_sorted)

#selected variables with optimal pi 
selected_vars1 <- selprop1[selprop1 >= 0.78]
print(selected_vars1)


#save selection proportions so SS doesnt have to be run again 
#selprop1_df <- data.frame(Variable = names(selprop1_sorted),SelectionProportion = selprop1_sorted)
#write.csv(selprop1_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selprop1_sorted.csv",row.names = FALSE)
#stable_vars_df <- data.frame(Variable = names(stable_vars_model1),SelectionProportion = stable_vars_model1)
#write.csv(stable_vars_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/stable_vars_model1.csv",row.names = FALSE)


# Visualisation of Stable Variables (pi >= 0.78)
#rename variables
selected_vars_df1 <- selected_vars_df1 %>%
  mutate(Variable = recode(Variable,
                           "sleep_pattern_Poor sleep pattern" = "Poor Sleep",
                           "iibs_2yr_Yes" = "IBS (2yr)",
                           "social_isolation_Isolated" = "Socially Isolated",
                           "edu_CSEs or equivalent" = "CSEs or Equivalent",
                           "employment_status_Unemployed" = "Unemployed",
                           "nucleated_red_blood_cell_count" = "Nucleated RBC Count",
                           "comp_body_10_Plumper" = "Plumper Body Type",
                           "HDL_cholesterol" = "HDL Cholesterol"
  ))
#bar plot
ggplot(selected_vars_df1, aes(x = reorder(Variable, SelectionProportion), y = SelectionProportion)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Selection Proportions of Stable Variables",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12)

#Visualisation of all variables (not just selection proportion >= 0.78)
selprop1_df$Variable <- gsub("_", " ", selprop1_df$Variable)

# Create and save plot
ggplot(selprop1_df, aes(x = reorder(Variable, SelectionProportion), y = SelectionProportion)) +
  geom_col(fill = "deeppink3") +
  geom_hline(yintercept = 0.78, linetype = "dashed", colour = "red", size = 0.5) +
  coord_flip() +
  labs(
    title = "Selection Proportions of All Variables",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 4)
  )

# Save as high-res PNG with enough height for 94 variables
ggsave(
  filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selection_plot_allvariables_model1.png",
  width = 10, height = 24, dpi = 300
)

################### LOGISTIC REGRESSION FOR MODEL 1 ############################
selected_vars <- c(
  "HDL_cholesterol",                 # biomarker
  "nucleated_red_blood_cell_count", # biomarker
  "comp_body_10",                   # categorical
  "social_isolation",               # categorical
  "iibs_2yr",                       # categorical
  "sleep_pattern",                  # categorical
  "employment_status",              # categorical
  "edu",                            # categorical
  "age",                            # numeric confounder
  "sex",                            # categorical confounder
  "ethnicity"                       # categorical confounder
)

# 1. Load datasets
X_train_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
X_train_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_exp_imputed_mice.rds")
Y_train     <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")
Y_train     <- as.numeric(as.character(Y_train))

# 2. Merge biomarker and raw exposure data
rownames(X_train_exp) <- rownames(X_train_bio)
X_train_exp$eid <- rownames(X_train_exp)
X_train_bio$eid <- rownames(X_train_bio)

X_full <- merge(X_train_bio, X_train_exp, by = "eid")
rownames(X_full) <- X_full$eid
X_full <- X_full[, !colnames(X_full) %in% "eid"]

# 3. Train-validation split
set.seed(2025)
n <- nrow(X_full)
train_size <- round(0.625 * n)
idx <- sample(1:n)
train_idx <- idx[1:train_size]
val_idx <- idx[(train_size + 1):n]

X_val <- X_full[val_idx, ]
Y_val <- Y_train[val_idx]

# 4. Build validation dataframe
df_val <- X_val[, selected_vars]
df_val$CKD <- Y_val

# 5. Relevel confounders for interpretability
df_val$sex <- relevel(as.factor(df_val$sex), ref = "Female")
df_val$ethnicity <- relevel(as.factor(df_val$ethnicity), ref = "White")

# 6. Fit logistic regression
fmla <- reformulate(selected_vars, response = "CKD")

#7. Calculate class weights
n_total <- nrow(df_val)
n_0 <- sum(df_val$CKD == 0)
n_1 <- sum(df_val$CKD == 1)
n_total <- nrow(df_val)

wts <- ifelse(df_val$CKD == 1, n_total / (2 * n_1), n_total / (2 * n_0))

logit_mod <- glm(fmla, data = df_val, family = binomial, weights = wts)
#logit_mod <- glm(fmla, data = df_val, family = binomial)

# 7. Output model summary
summary(logit_mod)

# 8. Odds ratios with 95% confidence intervals
se <- sqrt(diag(vcov(logit_mod)))
ORs <- exp(coef(logit_mod))
lower <- exp(coef(logit_mod) - 1.96 * se)
upper <- exp(coef(logit_mod) + 1.96 * se)

exp_coef <- data.frame(
  OR = round(ORs, 3),
  CI_lower = round(lower, 3),
  CI_upper = round(upper, 3)
)

# 9.  Forest Plot of Coefficients (ORs with CIs)
or_table <- as.data.frame(exp_coef)
or_table$Variable <- rownames(or_table)
colnames(or_table)[2:3] <- c("CI_low", "CI_high")

# Remove intercept for plotting and recode variables to make it pretty
or_table_plot <- or_table %>%
  filter(Variable != "(Intercept)") %>%
  mutate(Variable = factor(Variable, levels = Variable[order(OR)]))

or_table_plot <- or_table_plot %>%
  mutate(Variable = recode(Variable,
                           "HDL_cholesterol" = "HDL Cholesterol",
                           "nucleated_red_blood_cell_count" = "Nucleated RBC Count",
                           "comp_body_10Thinner" = "Thinner Body Type",
                           "comp_body_10Plumper" = "Plumper Body Type",
                           "social_isolationIsolated" = "Socially Isolated",
                           "iibs_2yrYes" = "IBS (2 yr)",
                           "sleep_patternIntermediate sleep pattern" = "Intermediate Sleep",
                           "sleep_patternPoor sleep pattern" = "Poor Sleep",
                           "employment_statusUnemployed" = "Employment: Unemployed",
                           "employment_statusRetired" = "Employment:Retired",
                           "employment_statusOther" = "Employment: Other Employment",
                           "employment_statusLooking after home and/or family" = "Employment: Looking After Home",
                           "employment_statusUnable to work due to sickness" = "Employment: Unable to Work",
                           "eduNone of the above" = "Education: No Qualifications",
                           "eduA levels/AS levels or equivalent" = "Education: A Levels",
                           "eduO levels/GCSEs or equivalent" = "Education: GCSEs",
                           "eduCSEs or equivalent" = "Education: CSEs",
                           "eduNVQ or HND or HNC or equivalent" = "Education: Vocational Qualification",
                           "eduOther professional qualifications eg: nursing, teaching" = "Education: Other Prof. Qualification",
                           "age" = "Age",
                           "sexMale" = "Male (vs Female)",
                           "ethnicityAsian" = "Asian",
                           "ethnicityBlack" = "Black",
                           "ethnicityOther" = "Other Ethnicity"
  ))


# Plot
ggplot(or_table_plot, aes(x = Variable, y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  coord_flip() +
  labs(title = "Odds Ratios with 95% CI (Model 1)",
       y = "Odds Ratio",
       x = "") +
  theme_minimal(base_size = 14)

saveRDS(or_table_plot,"/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/OR_model1.rds")
######### GET THE TOTAL DIRECT EFFECTS ############################
selected_exposures <- c(
  "sleep_pattern",
  "employment_status",
  "edu",
  "comp_body_10",
  "social_isolation",
  "iibs_2yr"
)
#Reference levels: "Healthy sleep pattern, Employed, College or University degree, average body, not isolated, no iibs in 2 years
# 1. Get the coefficient vector
coefs <- coef(logit_mod)

# 2. Extract all terms that match the selected exposures
direct_terms <- unlist(sapply(selected_exposures, function(var) {
  grep(paste0("^", var), names(coefs), value = TRUE)
}))

# 3. Sum up the direct effect terms (excluding reference levels)
total_direct_effect <- sum(coefs[direct_terms])
cat("Total Direct Effect:", round(total_direct_effect, 4), "\n") #Total Direct Effect: -0.0954

#4. Visualise
exposure_contributions <- lapply(selected_exposures, function(var) {
  terms <- grep(paste0("^", var), names(coefs), value = TRUE)
  total <- sum(coefs[terms])
  data.frame(Exposure = var, Total_Effect = total)
})
contrib_df <- do.call(rbind, exposure_contributions)
contrib_df <- contrib_df %>% arrange(desc(abs(Total_Effect)))

#write.csv(contrib_df,"/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/contribution_per_exposure_model1.csv",row.names = FALSE)

pretty_labels <- c(
  "sleep_pattern"      = "Sleep Pattern",
  "employment_status"  = "Employment Status",
  "edu"                = "Education Level",
  "comp_body_10"       = "Body Composition",
  "social_isolation"   = "Social Isolation",
  "iibs_2yr"           = "IBS (2-Year)"
)
contrib_df$Exposure <- recode(contrib_df$Exposure, !!!pretty_labels)

p <- ggplot(contrib_df, aes(x = reorder(Exposure, Total_Effect), y = Total_Effect)) +
  geom_col(fill = "navy") +
  coord_flip() +
  labs(
    title = "Total Direct Effect Contribution by Exposure",
    x = "Exposure Variable",
    y = "Sum of Coefficients"
  ) +
  theme_minimal()

# Save updated plot
ggsave(
  filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/contribution_plot_model1_named.png",
  plot = p,
  width = 8, height = 6, dpi = 300
)

