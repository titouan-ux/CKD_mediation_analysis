rm(list=ls())
project_path=dirname(rstudioapi::getActiveDocumentContext()$path)
##################
#dev.off()
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
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)


OR_direct_effect<- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/OR_model1.rds")
OR_direct_effect$Var <- rownames(OR_direct_effect)
OR_direct_effect <- OR_direct_effect %>% 
  mutate(logodd=log(OR),
         logodd_CI_low=log(CI_low),
         logodd_CI_high=log(CI_high),
         se= (log(CI_high) - log(CI_low)) / (2 * 1.96)
)
linear_RBC <- read_csv('/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/lm_coef_RBC_weights.csv',show_col_types = FALSE)
linear_RBC <- linear_RBC %>% mutate( se= (CI_high - CI_low) / (2 * 1.96))
linear_HDL <- read_csv('/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/lm_coef_HDL_weights.csv',show_col_types = FALSE)
linear_HDL <- linear_HDL %>% mutate( se= (CI_high - CI_low) / (2 * 1.96))

linear_RBC <- linear_RBC %>% mutate( indirect_effect=OR_direct_effect['nucleated_red_blood_cell_count','logodd'] * coef)
linear_HDL <- linear_HDL %>% mutate( indirect_effect=OR_direct_effect['HDL_cholesterol','logodd'] * coef)

total_indirect <- full_join(
  linear_RBC %>% select(Variable, indirect_effect) %>% rename(indirect_RBC = indirect_effect),
  linear_HDL %>% select(Variable, indirect_effect) %>% rename(indirect_HDL = indirect_effect),
  by = "Variable"
)
#Remove Intercept
total_indirect <- total_indirect[-c(1),]
total_indirect <- total_indirect %>%
  mutate(
#    indirect_RBC = replace_na(indirect_RBC, 0),
#    indirect_HDL = replace_na(indirect_HDL, 0),
    total_indirect_effect = ifelse(is.na(indirect_RBC) & is.na(indirect_HDL), NA, 
                                   ifelse(is.na(indirect_RBC), indirect_HDL, 
                                          ifelse(is.na(indirect_HDL), indirect_RBC, indirect_RBC+indirect_HDL))))
# Get logodd + se for biomarkers
logodd_rbc <- OR_direct_effect['nucleated_red_blood_cell_count', 'logodd']
logodd_rbc_se <- OR_direct_effect['nucleated_red_blood_cell_count', 'se']

logodd_hdl <- OR_direct_effect['HDL_cholesterol', 'logodd']
logodd_hdl_se <- OR_direct_effect['HDL_cholesterol', 'se']

# Add SE for each indirect effect
linear_RBC <- linear_RBC %>%
  mutate(
    indirect_se = sqrt((coef^2 * logodd_rbc_se^2) + (logodd_rbc^2 * se^2))
  )

linear_HDL <- linear_HDL %>%
  mutate(
    indirect_se = sqrt((coef^2 * logodd_hdl_se^2) + (logodd_hdl^2 * se^2))
  )

# Join SEs
total_indirect <- total_indirect %>%
  left_join(
    linear_RBC %>% select(Variable, se_RBC = indirect_se),
    by = "Variable"
  ) %>%
  left_join(
    linear_HDL %>% select(Variable, se_HDL = indirect_se),
    by = "Variable"
  ) %>%
  mutate(
#    se_RBC = replace_na(se_RBC, 0),
#    se_HDL = replace_na(se_HDL, 0),
    total_se = sqrt(se_RBC^2 + se_HDL^2),
    
    total_se = ifelse(is.na(se_RBC) & is.na(se_HDL), NA, 
                      ifelse(is.na(se_RBC), se_HDL, 
                             ifelse(is.na(se_HDL), se_RBC, sqrt(se_RBC^2 + se_HDL^2)))),
    CI_low_indirect = ifelse(is.na(total_indirect_effect),NA,total_indirect_effect - 1.96 * total_se),
    CI_high_indirect = ifelse(is.na(total_indirect_effect),NA,total_indirect_effect + 1.96 * total_se)
  )

total_indirect <- total_indirect %>%
  full_join(
    OR_direct_effect %>% select(Var, direct_effect = logodd, CI_low_direct=logodd_CI_low,CI_high_direct=logodd_CI_high),
    by = join_by("Variable"=="Var")
  )
#%>%
#  mutate(
#    direct_effect = replace_na(direct_effect, 0),
#    CI_low_direct = replace_na(CI_low_direct, 0),
#    CI_high_direct = replace_na(CI_high_direct, 0)
#    )
  
# Add CIs for individual indirect effects
linear_RBC <- linear_RBC %>%
  mutate(
    CI_low_indirect_RBC = indirect_effect - 1.96 * indirect_se,
    CI_high_indirect_RBC = indirect_effect + 1.96 * indirect_se
  )

linear_HDL <- linear_HDL %>%
  mutate(
    CI_low_indirect_HDL = indirect_effect - 1.96 * indirect_se,
    CI_high_indirect_HDL = indirect_effect + 1.96 * indirect_se
  )

# Now rejoin with full CIs for each path
total_indirect <- total_indirect %>%
  left_join(
    linear_RBC %>% select(Variable, CI_low_indirect_RBC, CI_high_indirect_RBC),
    by = "Variable"
  ) %>%
  left_join(
    linear_HDL %>% select(Variable, CI_low_indirect_HDL, CI_high_indirect_HDL),
    by = "Variable"
  ) %>%
  relocate(
    indirect_RBC, CI_low_indirect_RBC, CI_high_indirect_RBC,
    indirect_HDL, CI_low_indirect_HDL, CI_high_indirect_HDL,
    total_indirect_effect, CI_low_indirect, CI_high_indirect,
    direct_effect, CI_low_direct, CI_high_direct,
    .after = Variable
  )  
#%>%
#mutate(
#  CI_low_indirect_RBC = replace_na(CI_low_indirect_RBC, 0),
#  CI_high_indirect_RBC = replace_na(CI_high_indirect_RBC, 0),
#  CI_high_indirect_HDL = replace_na(CI_high_indirect_HDL, 0),
#  CI_low_indirect_HDL = replace_na(CI_low_indirect_HDL, 0)
#  )

#CHANGE NAMES
#total_indirect <- total_indirect %>% mutate(var_names = ifelse(Variable %in% rownames(OR_direct_effect),as.character(OR_direct_effect[OR_direct_effect$Var == Variable, "Variable"]),Variable))
total_indirect <- total_indirect %>%
  rowwise() %>%
  mutate(Variable = if (Variable %in% OR_direct_effect$Var) {
    as.character(OR_direct_effect$Variable[OR_direct_effect$Var == Variable])
  } else {
    Variable
  }) %>%
  ungroup()


######PLOT######
pretty_names <- c(
  "waist_circumference" = "Waist Circumference",
  "dist_greenspace" = "Distance to Greenspace",
  "own_rentOther" = "Housing: Other",
  "own_rentRent" = "Housing: Rented",
  "alcohol_categoryNever/Abstinent (Excluding previous drinkers)" = "Alcohol: Never/Abstinent",
  "alcohol_categoryModerate" = "Alcohol: Moderate",
  "alcohol_categoryHeavy" = "Alcohol: Heavy",
  "alcohol_categoryPrevious Drinker (who currently never drink)" = "Alcohol: Previous Drinker",
  "alcohol_categoryAbusive" = "Alcohol: Abusive",
  "ozone" = "Ozone",
  "PM10" = "PM10",
  "Age" = "Age",
  "Male (vs Female)" = "Sex: Male",
  "Asian" = "Ethnicity: Asian",
  "Black" = "Ethnicity: Black",
  "Other Ethnicity" = "Ethnicity: Other",
  "bmi" = "BMI",
  "sys_bp" = "Systolic Blood Pressure",
  "diet_score" = "Diet Score",
  "NO2" = "NO2",
  "imd_quin" = "Multiple Deprivation",
  "smoke_statusPrevious" = "Smoking: Previous",
  "smoke_statusCurrent" = "Smoking: Current",
  "Thinner Body Type" = "Body Type: Thinner",
  "Plumper Body Type" = "Body Type: Plumper",
  "comp_height_10Shorter" = "Height: Shorter",
  "comp_height_10Taller" = "Height: Taller",
  "physician_natdYes" = "Physician Diagnosed Depression",
  "Employment: Looking After Home" = "Employment: Home Duties",
  "Employment: Other Employment" = "Employment: Other",
  "Employment:Retired" = "Employment: Retired",
  "Employment: Unable to Work" = "Employment: Unable to Work",
  "Employment: Unemployed" = "Employment: Unemployed",
  "Education: No Qualifications" = "Education: None",
  "Education: A Levels" = "Education: A-Levels",
  "Education: GCSEs" = "Education: GCSEs",
  "Education: CSEs" = "Education: CSEs",
  "Education: Vocational Qualification" = "Education: Vocational",
  "Education: Other Prof. Qualification" = "Education: Other Qualification",
  "met_score" = "MET Score",
  "iibs_2yr"           = "IBS (2-Year)",
  "HDL_cholesterol" = "HDL Cholesterol",
  "nucleated_red_blood_cell_count" = "Nucleated RBC Count",
  "social_isolation"   = "Socially Isolated"
  
)

# Build a tidy data frame manually
plot_data <- bind_rows(
  total_indirect %>%
    select(Variable, value = total_indirect_effect, CI_low = CI_low_indirect, CI_high = CI_high_indirect) %>%
    mutate(effect_type = "Total Indirect"),
  total_indirect %>%
    select(Variable, value = direct_effect, CI_low = CI_low_direct, CI_high = CI_high_direct) %>%
    mutate(effect_type = "Direct")
) %>%
  # Remove rows where value == 0 and both CIs are 0 or NA
  filter(!(is.na(value) & (CI_low == 0 | is.na(CI_low)) & (CI_high == 0 | is.na(CI_high))))

plot_data <- plot_data %>%
  mutate(pretty_var = recode(Variable, !!!pretty_names))

# Step 1: Identify which variables have direct, indirect, or both
effect_presence <- total_indirect %>%
  mutate(
    has_direct = !is.na(direct_effect),
    has_indirect = !is.na(total_indirect_effect),
    effect_group = case_when(
      has_direct & !has_indirect ~ "Only Direct",
      !has_direct & has_indirect ~ "Only Indirect",
      has_direct & has_indirect ~ "Both Effects",
      TRUE ~ "Confounder" # this will catch any leftovers
    )
  ) %>%
  select(Variable, effect_group)
effect_presence[effect_presence$Variable %in% c("Age", "Male (vs Female)", "Asian", "Black", "Other Ethnicity"), "effect_group"]<-"Confounder"
effect_presence[effect_presence$Variable %in% c("Nucleated RBC Count","HDL Cholesterol"), "effect_group"]<-"Biomarkers"


plot_data <- plot_data %>%
  left_join(effect_presence, by = "Variable")


plot_data <- plot_data %>%
  group_by(effect_group) %>%
  mutate(var_order = dense_rank(desc(abs(value)))) %>%  # order by effect size within group
  ungroup() %>%
  arrange(factor(effect_group, levels = c("Only Indirect",  "Both Effects","Confounder", "Only Direct", "Biomarkers")), var_order)

custom_order <- c(
  "Education: None", "Education: GCSEs", "Education: CSEs", 
  "Education: A-Levels", "Education: Vocational", "Education: Other Qualification",
  "Employment: Home Duties", "Employment: Other", "Employment: Retired", "Employment: Unable to Work", "Employment: Unemployed",
  "Alcohol: Never/Abstinent", "Alcohol: Moderate", "Alcohol: Heavy", "Alcohol: Previous Drinker", "Alcohol: Abusive",
  "Smoking: Previous", "Smoking: Current",
  "Body Type: Thinner", "Body Type: Plumper",
  "Height: Shorter", "Height: Taller",
  "BMI", 
  "Waist Circumference", "Systolic Blood Pressure", "HDL Cholesterol", "Nucleated RBC Count",
  "Distance to Greenspace", "Housing: Rented", "Housing: Other",
  "NO2", "PM10", "Ozone",
  "Diet Score", "MET Score","Physician Diagnosed Depression","Multiple Deprivation",
  "Socially Isolated", "Poor Sleep", "Intermediate Sleep", "IBS (2 yr)",
  "Ethnicity: Asian", "Ethnicity: Black", "Ethnicity: Other","Age", "Sex: Male"
)

# Step 4: Make sure pretty_var is a factor with the desired order for plotting
plot_data$pretty_var <- factor(plot_data$pretty_var, levels = custom_order)

ggplot(plot_data[plot_data$effect_group %in% c("Both Effects", "Confounder"),], 
       aes(x = pretty_var , y = exp(value), color = effect_type)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = exp(CI_low), ymax = exp(CI_high)), width = 0.2,
                position = position_dodge(width = 0.6)) +
  coord_flip() +
  theme_minimal(base_size = 18) +  # Increase base size for all text
  theme(
    axis.text = element_text(size = 18),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
    legend.title = element_text(size = 15),    # Legend title
    legend.text = element_text(size = 14),     # Legend text
    plot.title = element_text(size = 18, face = "bold")  # Plot title
  ) +
  labs(x = NULL, y = "Odds ratios",
       title = NULL,
       color = "Effect Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Direct" = "red",
    "Total Indirect" = "blue"
  ))

ggplot(plot_data[plot_data$effect_group %in% c("Only Direct", "Biomarkers"),], 
       aes(x = pretty_var , y = exp(value), color = effect_type), xlim) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = exp(CI_low), ymax = exp(CI_high)), width = 0.2,
                position = position_dodge(width = 0.6)) +
  coord_flip() +
  theme_minimal(base_size = 18) +  # Increase base size for all text
  theme(
    axis.text = element_text(size = 18),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
    legend.title = element_text(size = 15),    # Legend title
    legend.text = element_text(size = 14),     # Legend text
    plot.title = element_text(size = 18, face = "bold")  # Plot title
  ) +
  labs(x = NULL, y = "Odds ratios",
       title = NULL,
       color = "Effect Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Direct" = "red",
    "Total Indirect" = "blue"
  ))

ggplot(plot_data[plot_data$effect_group %in% c("Only Indirect"),], 
       aes(x = pretty_var , y = exp(value), color = effect_type), xlim) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = exp(CI_low), ymax = exp(CI_high)), width = 0.2,
                position = position_dodge(width = 0.6)) +
  coord_flip() +
  theme_minimal(base_size = 18) +  # Increase base size for all text
  theme(
    axis.text = element_text(size = 18),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
    legend.title = element_text(size = 15),    # Legend title
    legend.text = element_text(size = 14),     # Legend text
    plot.title = element_text(size = 18, face = "bold")  # Plot title
  ) +
  labs(x = NULL, y = "Odds ratios",
       title = NULL,
       color = "Effect Type") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "Direct" = "red",
    "Total Indirect" = "blue"
  ))


############HEATMAP###############
heatmap_data <- bind_rows(
  total_indirect %>%
    select(Variable, value = indirect_RBC) %>%
    mutate(effect_type = "Indirect via RBC"),
  total_indirect %>%
    select(Variable, value = indirect_HDL) %>%
    mutate(effect_type = "Indirect via HDL"),
  total_indirect %>%
    select(Variable, value = total_indirect_effect) %>%
    mutate(effect_type = "Total Indirect"),
  total_indirect %>%
    select(Variable, value = direct_effect) %>%
    mutate(effect_type = "Direct")
)

heatmap_data <- heatmap_data %>%
  mutate(pretty_var = recode(Variable, !!!pretty_names))
heatmap_data$pretty_var <- factor(heatmap_data$pretty_var, levels = custom_order)

heatmap_data <- heatmap_data[-which(heatmap_data$Variable %in% c("Nucleated RBC Count","HDL Cholesterol")),]

max_abs_val <- max(abs(heatmap_data$value), na.rm = TRUE)
# Step 3: Create heatmap
ggplot(heatmap_data, aes(x = effect_type, y = pretty_var, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", direction = 1,limits = c(-max_abs_val, max_abs_val),
                       oob = scales::squish
  )+
#scale_fill_viridis_c(option = "viridis", na.value = "grey80", direction = 1) +  
  labs(x = NULL, y = NULL, fill = "Effect\n(log-odds)",
       title =NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),       # Axis tick labels
    legend.title = element_text(size = 15),    # Legend title
    legend.text = element_text(size = 14),     # Legend text
    plot.title = element_text(size = 14, face = "bold"),  # Plot title
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Classify variables
heatmap_data <- heatmap_data %>%
  mutate(var_type = ifelse(grepl(":", pretty_var), "Categorical", "Continuous"))

heatmap_data[which(heatmap_data$Variable=="Socially Isolated"),"var_type"] = "Categorical"
heatmap_data[which(heatmap_data$Variable=="Poor Sleep"),"var_type"] = "Categorical"
heatmap_data[which(heatmap_data$Variable=="Intermediate Sleep"),"var_type"] = "Categorical"
# Now split the data
heatmap_data_cat <- heatmap_data %>% filter(var_type == "Categorical")
heatmap_data_cont <- heatmap_data %>% filter(var_type == "Continuous")

max_abs_val <- max(abs(heatmap_data_cont$value), na.rm = TRUE)
ggplot(heatmap_data_cont, aes(x = effect_type, y = pretty_var, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", direction = 1,limits = c(-max_abs_val, max_abs_val),
                       oob = scales::squish
  )+
  #scale_fill_viridis_c(option = "viridis", na.value = "grey80", direction = 1) +  
  labs(x = NULL, y = NULL, fill = "Effect\n(log-odds)",
       title =NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),       # Axis tick labels
    legend.title = element_text(size = 15),    # Legend title
    legend.text = element_text(size = 14),     # Legend text
    plot.title = element_text(size = 14, face = "bold"),  # Plot title
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
#####BARPPLOT
# Prepare the full data
barplot_data <- bind_rows(
  total_indirect %>%
    select(Variable, value = indirect_RBC) %>%
    mutate(effect_type = "Indirect via RBC"),
  total_indirect %>%
    select(Variable, value = indirect_HDL) %>%
    mutate(effect_type = "Indirect via HDL"),
  total_indirect %>%
    select(Variable, value = total_indirect_effect) %>%
    mutate(effect_type = "Total Indirect"),
  total_indirect %>%
    select(Variable, value = direct_effect) %>%
    mutate(effect_type = "Direct")
)

barplot_data <- barplot_data %>%
  mutate(pretty_var = recode(Variable, !!!pretty_names))
barplot_data$pretty_var <- factor(barplot_data$pretty_var, levels = custom_order)

# Remove unwanted variables
barplot_data <- barplot_data[-which(barplot_data$Variable %in% c("Nucleated RBC Count", "HDL Cholesterol")),]

# Classify variables
barplot_data <- barplot_data %>%
  mutate(var_type = ifelse(grepl(":", pretty_var), "Categorical", "Continuous"))

barplot_data[which(barplot_data$Variable=="Socially Isolated"),"var_type"] = "Categorical"
barplot_data[which(barplot_data$Variable=="Poor Sleep"),"var_type"] = "Categorical"
barplot_data[which(barplot_data$Variable=="Intermediate Sleep"),"var_type"] = "Categorical"

# Now split the data
barplot_data_cat <- barplot_data %>% filter(var_type == "Categorical")
barplot_data_cont <- barplot_data %>% filter(var_type == "Continuous")

# Define color scheme
effect_colors <- c(
  "Direct" = "#1f77b4",            # Blue
  "Total Indirect" = "#ff7f0e",    # Orange
  "Indirect via HDL" = "#2ca02c",  # Green
  "Indirect via RBC" = "#d62728"   # Red
)

# Plot for Categorical Variables
p_cat <- ggplot(barplot_data_cat, aes(x = pretty_var, y = value, fill = effect_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = effect_colors) +
  labs(x = NULL, y = "Effect (log-odds)", fill = "Effect Type",
       title = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  coord_flip()


# Plot for Continuous Variables
p_cont <- ggplot(barplot_data_cont, aes(x = pretty_var, y = value, fill = effect_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = effect_colors) +
  labs(x = NULL, y = "Effect (log-odds)", fill = "Effect Type",
       title = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  coord_flip()

# Display
p_cat
p_cont


library(tidyr)
library(dplyr)

barplot_table <- barplot_data %>%
  pivot_wider(names_from = effect_type, values_from = value) %>%
  mutate(
    Total = Direct + `Total Indirect`,
    `Direct Effect (%)` = 100 * Direct / Total,
    `Indirect Effect (%)` = 100 * `Total Indirect` / Total
  ) %>%
  select(Variable, `Direct Effect (%)`, `Indirect Effect (%)`)

barplot_table <- barplot_data %>%
  pivot_wider(names_from = effect_type, values_from = value) %>%
  mutate(
    Direct = ifelse(is.na(Direct), 0, Direct),
    `Total Indirect` = ifelse(is.na(`Total Indirect`), 0, `Total Indirect`)
  ) %>%
  mutate(
    Total = abs(Direct) + abs(`Total Indirect`),
    `Direct Effect (%)` = 100 * abs(Direct) / Total,
    `Indirect Effect (%)` = 100 * abs(`Total Indirect`) / Total
  ) %>%
  select(Variable, `Direct Effect (%)`, `Indirect Effect (%)`)


