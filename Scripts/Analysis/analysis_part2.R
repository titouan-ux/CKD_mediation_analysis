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
# library(missForestPredict)
library(caret)
library(pROC)
library(fake)
library(igraph)
library(pheatmap)
library(sharp)
# library(onehot)

#load data 

selected_vars <- read_csv("/rds/general/user/yc8024/projects/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars1.csv")
selected_vars <- selected_vars[[1]]
selected_bio <- c("HDL_cholesterol", "nucleated_red_blood_cell_count")
# selected_exp <- selected_vars[!selected_vars %in% selected_bio]


X_train_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
X_test_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_biomarker.rds")
X_train_bio<- X_train_bio[, selected_bio]
X_test_bio<- X_test_bio[, selected_bio]

X_train_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_exposure.rds")
X_test_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_test_exposure.rds")
# X_train_exp <- X_train_exp[, selected_exp]
# X_test_exp <- X_test_exp[, selected_exp]

#nb y vectors are the same for exposures and biomarkers

# Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")
# Y_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_test.rds")

# remove age and sex and ethnicity from exposures to prevent leakage, these will be added back as confounders later 
## DO WE WANT TO REMOVE ANY OTHERS AS CONFOUNDERS?? 

#X_train_exp <- X_train_exp[, !(colnames(X_train_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]
#X_test_exp  <- X_test_exp[, !(colnames(X_test_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]

# add eid to biomarker columns 

# X_train_bio$eid <- rownames(X_train_bio)
# X_test_bio$eid <- rownames(X_test_bio)

# Merge datasets by 'eid' column
# X_train_all <- merge(X_train_bio, X_train_exp, by = "eid")
# X_train_all <- X_train_all[, !colnames(X_train_all) %in% c("eid")]
# 
# X_test_all <- merge(X_test_bio, X_test_exp, by = "eid")
# X_test_all <- X_test_all[, !colnames(X_test_all) %in% c("eid")]


## split train_all into train_all 50% and validation 30% of the dataset 
set.seed(2025)

n <- nrow(X_train_exp)
train_prop <- 0.625
val_prop <- 0.375
train_size <- round(train_prop * n)

# Step 3: Sample indices
all_indices <- sample(1:n)
train_indices <- all_indices[1:train_size]
val_indices <- all_indices[(train_size + 1):n]

# Step 4: Create data splits
#X_val_exp   <- X_train_exp[val_indices, ]
X_train_exp <- X_train_exp[train_indices, ]

#Y_val_bio <- X_train_bio[val_indices, ]
Y_train_bio <- X_train_bio[train_indices, ]

###################################

# format data for analysis 

# Convert X data frames to matrices

X_train_exp <- as.matrix(X_train_exp)
X_test_exp <- as.matrix(X_test_exp)

Y_train_bio1 <- as.numeric(Y_train_bio[[1]])
Y_test_bio1 <- as.matrix(X_test_bio[[1]])

Y_train_bio2 <- as.numeric(Y_train_bio[[2]])
Y_test_bio2 <- as.matrix(X_test_bio[[2]])

# X_train_all <- as.matrix(X_train_all)
# X_test_all <- as.matrix(X_test_all)

# X_val_all <- as.matrix(X_val_all)

# Convert Y to numeric vectors
# Y_train <- as.numeric(as.character(Y_train))
# Y_test <- as.numeric(as.character(Y_test))

#####################

#LASSO stability selection for indirect model (i.e. using all exposures and outcome joint matrix)
## n.b. do we want to constrain calibration with PFER
# ---------------------- biomarker = HDL_cholesterol
set.seed(2025)

X_conf <- X_train_exp[, colnames(X_train_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black")]
X_train_exp  <- X_train_exp[, !(colnames(X_train_exp) %in% c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black"))]

t0 <- Sys.time()

model1 <- VariableSelection(
  xdata = cbind(X_train_exp, X_conf),
  ydata = Y_train_bio1,
  family = "gaussian",
  n_cat = 3,
  pi_list = seq(0.5, 0.99, by = 0.01), 
  penalty.factor = c(rep(1,
                         ncol(X_train_exp)), rep(0, ncol(X_conf)))
)
# model1 <- VariableSelection(
#   xdata = X_train_exp,
#   ydata = Y_train_bio1,
#   family = "gaussian",
#   n_cat = 3,
#   pi_list = seq(0.5, 0.99, by = 0.01)
# )
t1 <- Sys.time()
print(t1 - t0)

CalibrationPlot(model1)
# png("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/calibration_plot_HDL.png",
#     width = 1200, height = 800, res = 150)
# par(mar = c(6, 6, 6, 6))  # optional: adjust margins if needed
# CalibrationPlot(model1)
# dev.off()

selprop_HDL <- SelectionProportions(model1)
selprop_HDL <- sort(selprop_HDL, decreasing = TRUE)
# print(selprop_HDL)

# Calibrated parameters
optimal_par1 <- Argmax(model1) # lambda = 0.008985489; pi = 0.99
print(optimal_par1)

# Visualisation of selection proportions
par(mar = c(10, 5, 1, 1))
plot(selprop_HDL, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop_HDL >=
                                                               optimal_par1[2], yes = "red", no = "grey"), cex.lab = 1.5)
abline(h = optimal_par1[2], lty = 2, col = "darkred")
for (i in 1:length(selprop_HDL)) {
  axis(side = 1, at = i, labels = names(selprop_HDL)[i],
       las = 2, col = ifelse(selprop_HDL[i] >= optimal_par1[2],
                             yes = "red", no = "grey"), col.axis = ifelse(selprop_HDL[i] >=
                                                                            optimal_par1[2], yes = "red", no = "grey"))
}

#see what vars are selected at the calibrated pi 
stable_vars1 <- selprop_HDL[selprop_HDL >= 0.99]
print(stable_vars1)
stable_vars1_df <- as.data.frame(stable_vars1)
stable_vars1_df$variable <- rownames(stable_vars1_df)
names(stable_vars1_df)[1] <- c("selection_proportion")

#save selection proportions so SS doesnt have to be run again 
selprop_HDL_df <- as.data.frame(selprop_HDL)
selprop_HDL_df$variable <- rownames(selprop_HDL_df)
names(selprop_HDL_df)[1] <- c("selection_proportion")

# write.csv(selprop_HDL_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/sel_prop_HDL.csv", row.names = FALSE)
# write.csv(stable_vars1_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_HDL.csv", row.names = FALSE)


# ----------------------------- biomarker = nucleated_red_blood_cell_count
t0 <- Sys.time()
model2 <- VariableSelection(
  xdata = X_train_exp,
  ydata = Y_train_bio2,
  family = "gaussian",
  n_cat = 3,
  pi_list = seq(0.5, 0.99, by = 0.01)
)

# model2 <- VariableSelection(
#   xdata = cbind(X_train_exp, X_conf),
#   ydata = Y_train_bio2,
#   family = "gaussian",
#   n_cat = 3,
#   pi_list = seq(0.5, 0.99, by = 0.01), 
#   penalty.factor = c(rep(1,
#                          ncol(X_train_exp)), rep(0, ncol(X_conf)))
# )

t1 <- Sys.time()
print(t1 - t0)


CalibrationPlot(model2)
# png("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/calibration_plot_RBC.png", 
#     width = 1200, height = 800, res = 150)
# par(mar = c(6, 6, 6, 6))  # optional: adjust margins if needed
# CalibrationPlot(model2)
# dev.off()

selprop_RBC <- SelectionProportions(model2)
selprop_RBC <- sort(selprop_RBC, decreasing = TRUE)
# print(selprop_RBC)

# Calibrated parameters
optimal_par2 <- Argmax(model2) # lambda = 0.006042341, pi = 0.97
print(optimal_par2)

# Visualisation of selection proportions
par(mar = c(10, 5, 1, 1))
plot(selprop_RBC, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop_RBC >=
                                                               optimal_par2[2], yes = "red", no = "grey"), cex.lab = 1.5)
abline(h = optimal_par2[2], lty = 2, col = "darkred")
for (i in 1:length(selprop_RBC)) {
  axis(side = 1, at = i, labels = names(selprop_RBC)[i],
       las = 2, col = ifelse(selprop_RBC[i] >= optimal_par1[2],
                             yes = "red", no = "grey"), col.axis = ifelse(selprop_RBC[i] >=
                                                                            optimal_par2[2], yes = "red", no = "grey"))
}

#see what vars are selected at the calibrated pi 
stable_vars2 <- selprop_RBC[selprop_RBC >= 0.97]
print(stable_vars2)
stable_vars2_df <- as.data.frame(stable_vars2)
stable_vars2_df$variable <- rownames(stable_vars2_df)
names(stable_vars2_df)[1] <- c("selection_proportion")

# save selection proportions so SS doesnt have to be run again 
selprop_RBC_df <- as.data.frame(selprop_RBC)
selprop_RBC_df$variable <- rownames(selprop_RBC_df)
names(selprop_RBC_df)[1] <- c("selection_proportion")

# write.csv(selprop_RBC_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/sel_prop_RBC.csv", row.names = FALSE)
# write.csv(stable_vars2_df,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_RBC.csv", row.names = FALSE)


## visulize stable selected variables ----------------------
## HDL --------------
stable_vars_HDL <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_HDL.csv")
vars_HDL <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/sel_prop_HDL.csv")

#rename variables
stable_vars_HDL <- stable_vars_HDL %>%
  mutate(variable = recode(variable,
                           "waist_circumference" = "Waist Circumference",
                           "bmi"                       = "BMI",
                           "sys_bp"                    = "Systolic Blood Pressure",
                           "diet_score"                = "Diet Score",
                           "NO2"                       = "NO2",
                           "ozone"                     = "Ozone",
                           "imd_quin"                  = "Multiple Deprivation",
                           "alcohol_category_Never/Abstinent (Excluding previous drinkers)" = "Alcohol: Never/Abstinent",
                           "alcohol_category_Moderate" = "Alcohol: Moderate",
                           "alcohol_category_Heavy" = "Alcohol: Heavy",
                           "alcohol_category_Previous Drinker (who currently never drink)" =  "Alcohol: Previous Drinker",
                           "alcohol_category_Abusive" = "Alcohol: Abusive",
                           "smoke_status_Current" = "Smoking: Current",
                           "comp_body_10_Thinner" = "Body Type: Thinner",
                           "comp_body_10_Plumper" = "Body Type: Plumper",
                           "comp_height_10_Shorter" = "Height: Shorter",
                           "comp_height_10_Taller" = "Height: Taller",
                           "physician_natd_Yes" = "Physician Diagnosed Depression",
                           "employment_status_Looking after home and/or family" = "Employment: Home Duties",
                           "employment_status_Retired" = "Employment: Retired",
                           "edu_NVQ or HND or HNC or equivalent" = "Education: Other Qualification",
                           "met_score" = "MET Score"
  ))


## rename all variables
vars_HDL <-  vars_HDL%>%
  mutate(variable = case_when(
    # core numeric vars
    variable == "waist_circumference" ~ "Waist Circumference",
    variable == "bmi"                  ~ "BMI",
    variable == "sys_bp"               ~ "Systolic Blood Pressure",
    variable == "diet_score"           ~ "Diet Score",
    variable == "NO2"                  ~ "NO2",
    variable == "ozone"                ~ "Ozone",
    variable == "imd_quin"             ~ "Multiple Deprivation",
    variable == "met_score"            ~ "MET Score",
    variable == "impervious"           ~ "Impervious Surface Cover",
    variable == "neuroticism_score"    ~ "Neuroticism Score",
    # alcohol categories
    str_starts(variable, "alcohol_category_") ~
      paste0("Alcohol: ", str_to_title(str_replace(variable, "alcohol_category_", "") %>% 
                                         str_replace_all("_", " "))),
    # smoking categories
    str_starts(variable, "smoke_status_") ~
      paste0("Smoking: ", str_to_title(str_replace(variable, "smoke_status_", "") %>% 
                                                str_replace_all("_", " "))),
    # body-composition / height percentiles
    str_starts(variable, "comp_body_10_") ~
      paste0("Body Type: ",
             str_to_title(str_replace(variable, "comp_body_10_", "") %>% 
                            str_replace_all("_", " "))),
    str_starts(variable, "comp_height_10_") ~
      paste0("Height: ",
             str_to_title(str_replace(variable, "comp_height_10_", "") %>% 
                            str_replace_all("_", " "))),
    # physician-noted diagnosis
    variable == "physician_natd_Yes" ~ "Physician Diagnosed Depression",
    # employment
    str_starts(variable, "employment_status_") ~
      paste0("Employment: ",
             str_to_title(str_replace(variable, "employment_status_", "") %>% 
                            str_replace_all("_", " "))),
    # education
    str_starts(variable, "edu_") ~
      paste0("Education: ",
             str_to_title(str_replace(variable, "edu_", "") %>% 
                            str_replace_all("_| or ", " "))),
    # housing
    variable %in% c("own_rent_Other", "own_rent_Rent") ~
      paste0("Housing: ", str_to_title(str_replace(variable, "own_rent_", ""))),
    # income strata
    str_starts(variable, "avg_hhincome_") ~
      paste0("Household Income: ",
             str_to_title(str_replace(variable, "avg_hhincome_", "") %>% 
                            str_replace_all("_to_|_", " to "))),
    # sleep patterns
    str_starts(variable, "sleep_pattern_") ~
      paste0("Sleep Pattern: ",
             str_to_title(str_replace(variable, "sleep_pattern_", "") %>% 
                            str_replace_all("_", " "))),
    # fallback: title-case any remaining snake_case
    TRUE ~ str_to_title(str_replace_all(variable, "_", " "))
  ))
  
# visualisation of selected vars
g_HDL <- ggplot(stable_vars_HDL, aes(x = reorder(variable, selection_proportion), y = selection_proportion)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Selection Proportions of Stable Variables (HDL)",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12)

# Visualisation of all variables (not just selection proportion >= 0.99)
# Create and save plot
g_vars_HDL <- ggplot(vars_HDL, aes(x = reorder(variable, selection_proportion), y = selection_proportion)) +
  geom_col(fill = "deeppink3") +
  geom_hline(yintercept = 0.99, linetype = "dashed", colour = "red", size = 0.5) +
  coord_flip() +
  labs(
    title = "Selection Proportions of All Variables",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12)
  )


# Save as high-res PNG with enough height for 94 variables
# ggsave(g_HDL,
#   filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selection_plot_stable_vars_HDL.png",
#   width = 10, height = 6, dpi = 300
# )
# 
# ggsave(g_vars_HDL,
#        filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selection_plot_all_vars_HDL.png",
#        width = 10, height = 24, dpi = 300
# )



## RBC ------------------------

stable_vars_RBC <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_RBC.csv")
vars_RBC <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/sel_prop_RBC.csv")

stable_vars_RBC <- stable_vars_RBC %>%
  mutate(variable = recode(variable,
                           "waist_circumference" = "Waist Circumference",
                           "dist_greenspace" = "Distance to Greenspace",
                           "own_rent_Rent" = "Housing: Rented",
                           "alcohol_category_Never/Abstinent (Excluding previous drinkers)" = "Alcohol: Never/Abstinent",
                           "alcohol_category_Heavy" = "Alcohol: Heavy",
                           "ozone"                     = "Ozone",
                           "PM10" = "PM10"
  ))

# visualisation of selected vars
g_RBC <- ggplot(stable_vars_RBC, aes(x = reorder(variable, selection_proportion), y = selection_proportion)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Selection Proportions of Stable Variables (RBC)",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12)


# Visualisation of all variables (not just selection proportion >= 0.99)
# Create and save plot
g_vars_RBC <- ggplot(vars_RBC, aes(x = reorder(variable, selection_proportion), y = selection_proportion)) +
  geom_col(fill = "deeppink3") +
  geom_hline(yintercept = 0.97, linetype = "dashed", colour = "red", size = 0.5) +
  coord_flip() +
  labs(
    title = "Selection Proportions of All Variables (RBC)",
    x = "Variables",
    y = "Selection Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12)
  )

g_vars_HDL

# Save as high-res PNG with enough height for 94 variables
ggsave(g_RBC,
       filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selection_plot_stable_vars_RBC.png",
       width = 10, height = 6, dpi = 300
)

ggsave(g_vars_RBC,
       filename = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selection_plot_all_vars_RBC.png",
       width = 10, height = 24, dpi = 300
)




# ------------------------------- linear regression, HDL
vars_HDL <- c("waist_circumference",
              "bmi",
              "sys_bp",
              "diet_score",
              "NO2",
              "ozone",
              "imd_quin",
              "alcohol_category",
              "smoke_status",
              "comp_body_10",
              "comp_height_10",
              "physician_natd",
              "employment_status",
              "edu",
              "met_score",
              "age",
              "sex",
              "ethnicity")

vars_RBC <- c("waist_circumference",
              "dist_greenspace",
              "own_rent",
              "alcohol_category",
              "ozone",
              "PM10",
              "age",
              "sex",
              "ethnicity")

# 1. Load datasets

X_exp <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/Imputed_data/X_train_exp_imputed_mice.rds")
X_train_bio <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
selected_bio <- c("HDL_cholesterol", "nucleated_red_blood_cell_count")
X_train_bio <- X_train_bio[, selected_bio]

Y_train     <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")


# 2. Train-validation split
set.seed(2025)
n <- nrow(X_exp)
train_prop <- 0.625
val_prop <- 0.375
train_size <- round(train_prop * n)

# sample indices
all_indices <- sample(1:n)
train_indices <- all_indices[1:train_size]
val_indices <- all_indices[(train_size + 1):n]

X_val_exp <- X_exp[val_indices, ]
Y_val_bio <- X_train_bio[val_indices, ]

Y_val_bio1 <- Y_val_bio[[1]]
Y_val_bio2 <- as.numeric(Y_val_bio[[2]])

# 3. Build validation dataframe
df_val_HDL <- X_val_exp[, vars_HDL]
df_val_HDL$HDL <- Y_val_bio1

df_val_RBC <- X_val_exp[, vars_RBC]
df_val_RBC$RBC <- Y_val_bio2

# 4. Relevel confounders for interpretability
# df_val_RBC$sex <- relevel(as.factor(df_val$sex), ref = "Female")
# df_val$ethnicity <- relevel(as.factor(df_val$ethnicity), ref = "White")

# Calculate class weights
Y_val_CKD <- Y_train[val_indices]

n_total <- nrow(Y_val_CKD)
n_0 <- sum(Y_val_CKD == 0)
n_1 <- sum(Y_val_CKD == 1)
n_total <- length(Y_val_CKD)
wts_HDL <- ifelse(Y_val_CKD == 1, n_total / (2 * n_1), n_total / (2 * n_0))


# 5. Fit linear regression
lm_HDL <- lm(HDL ~ ., data = df_val_HDL, weights = wts_HDL)
summary(lm_HDL)

se <- sqrt(diag(vcov(lm_HDL)))
coefs <- coef(lm_HDL)
lower <- coef(lm_HDL) - 1.96 * se
upper <- coef(lm_HDL) + 1.96 * se

coef_table <- data.frame(
  coef = round(coefs, 6),
  CI_lower = round(lower, 6),
  CI_upper = round(upper, 6)
)

coef_table_HDL <- as.data.frame(coef_table)
coef_table_HDL$Variable <- rownames(coef_table_HDL)
colnames(coef_table_HDL)[2:3] <- c("CI_low", "CI_high")


# ---------------
lm_RBC <- lm(RBC ~ ., data = df_val_RBC, weights = wts_HDL)
summary(lm_RBC)

se <- sqrt(diag(vcov(lm_RBC)))
coefs <- coef(lm_RBC)
lower <- coef(lm_RBC) - 1.96 * se
upper <- coef(lm_RBC) + 1.96 * se

coef_table <- data.frame(
  coef = round(coefs, 6),
  CI_lower = round(lower, 6),
  CI_upper = round(upper, 6)
)

coef_table_RBC <- as.data.frame(coef_table)
coef_table_RBC$Variable <- rownames(coef_table_RBC)
colnames(coef_table_RBC)[2:3] <- c("CI_low", "CI_high")


write.csv(coef_table_HDL,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/lm_coef_HDL_weights.csv", row.names = FALSE)
write.csv(coef_table_RBC,file = "/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/lm_coef_RBC_weights.csv", row.names = FALSE)

###------------------------- Step 9---------------------------------

# Load all selected variables
selected_vars1 <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars1.csv")
selected_vars_HDL <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_HDL.csv")
selected_vars_RBC <- read.csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_RBC.csv")

# Load datasets
Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")
Y_test <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_test.rds")

# Combine all selected exposures
all_selected_exp <- c(as.character(selected_vars1$variable),
                      as.character(selected_vars_HDL$variable),
                      as.character(selected_vars_RBC$variable),
                      "social_isolation_Isolated", "sleep_pattern_Intermediate sleep pattern", "sleep_pattern_Poor sleep pattern", "iibs_2yr_Yes"
)
all_selected_exp <- unique(all_selected_exp)
list(all_selected_exp)

# Add confounders
confounders <- c("age", "sex_Male", "ethnicity_Asian", "ethnicity_Other", "ethnicity_Black")
selected_biomarkers <- c("HDL_cholesterol", "nucleated_red_blood_cell_count")

# train/validation split for the origin datasets
set.seed(2025)

# Load original dataset
X_train_exposure <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_exposure.rds")
X_train_biomarker <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/X_train_biomarker.rds")
Y_train <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/analysis_data_imputed_scaled/Y_train.rds")

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

# validation sets
X_val30_exp <- X_train_exposure[val_indices, ]
X_val30_bio <- X_train_biomarker[val_indices, ]
Y_val30 <- Y_train[val_indices]

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

