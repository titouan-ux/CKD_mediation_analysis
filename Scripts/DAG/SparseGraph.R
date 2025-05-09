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
################### SPARSE GRAPH ########################### 

X_exposures <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_external.rds"))
exposures <- colnames(X_exposures)
X_biomarkers <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_internal.rds"))
biomarkers <- colnames(X_biomarkers)
confounders <- c("age", "sex", "ethnicity") 
outcome <- c("CKD")

OR_model1 <- readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/OR_model1.rds")

selected_bio <- rownames(OR_model1)[rownames(OR_model1) %in% biomarkers]

rownames(OR_model1)[!rownames(OR_model1) %in% biomarkers]

selected_exposures_outcome <- c(
  "sleep_pattern",
  "employment_status",
  "edu",
  "comp_body_10",
  "social_isolation",
  "iibs_2yr"
)

coef_model2_HDL <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_HDL.csv",show_col_types = FALSE)
coef_model2_HDL$variable

selected_exposures_bio_HDL <- c("waist_circumference", "sys_bp","NO2","imd_quin","alcohol_category","smoke_status", "comp_body_10", "comp_height_10",
                                "employment_status","edu","bmi","ozone","physician_natd","met_score")

coef_model2_RBC <- read_csv("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/selected_vars_RBC.csv",show_col_types = FALSE)
coef_model2_RBC$variable

selected_exposures_bio_RBC <- c("waist_circumference", "own_rent","alcohol_category", "dist_greenspace","PM10", "ozone")

only_direct_exposures = selected_exposures_outcome[ ! selected_exposures_outcome %in% union(selected_exposures_bio_HDL,selected_exposures_bio_RBC)]
all <- intersect(intersect(selected_exposures_bio_HDL, selected_exposures_outcome),selected_exposures_bio_RBC) 
only_RBC <- selected_exposures_bio_RBC[!selected_exposures_bio_RBC %in% union(selected_exposures_bio_HDL, selected_exposures_outcome)]
only_HDL <- selected_exposures_bio_HDL[!selected_exposures_bio_HDL %in% union(selected_exposures_bio_RBC, selected_exposures_outcome)]
RBD_HDL <- intersect(selected_exposures_bio_HDL,selected_exposures_bio_RBC)[ ! intersect(selected_exposures_bio_HDL,selected_exposures_bio_RBC) %in% selected_exposures_outcome]
RBC_outcome <- intersect(selected_exposures_outcome,selected_exposures_bio_RBC)[ ! intersect(selected_exposures_outcome,selected_exposures_bio_RBC) %in% selected_exposures_bio_HDL]
HDL_outcome <- intersect(selected_exposures_outcome,selected_exposures_bio_HDL)[ ! intersect(selected_exposures_outcome,selected_exposures_bio_HDL) %in% selected_exposures_bio_RBC]
exposures <- c(only_direct_exposures,HDL_outcome,only_HDL,all,RBC_outcome,RBD_HDL,only_RBC)


# Create a list of all nodes
all_nodes <- c(confounders, selected_bio, outcome, exposures)

node_category <- c(
  rep("confounder", length(confounders)),
  rep("biomarker", length(selected_bio)),
  rep("outcome", length(outcome)),
  rep("direct_exposure", length(only_direct_exposures)),
  rep("mixed_exposure", length(c(HDL_outcome,RBC_outcome, all,RBD_HDL))),
  rep("indirect_exposure", length(c(only_HDL,only_RBC)))
)
# Define edges for the DAG - only between different categories
edges <- c()
# Create the graph with vertices but no edges
g <- make_empty_graph(n = length(all_nodes), directed = TRUE)
# Set node attributes
V(g)$name <- all_nodes
V(g)$category <- node_category

for (b in selected_bio) {
  g <- add_edges(g, c(b,outcome))
}
for (e in selected_exposures_outcome) {
  g <- add_edges(g, c(e, outcome))
}
for (c in confounders) {
  g <- add_edges(g, c(c, outcome))
}
for (c in confounders) {
  for (b in selected_bio){
    g <- add_edges(g, c(c, b))
  }
}
for (e in selected_exposures_bio_HDL) {
  g <- add_edges(g, c(e, "HDL_cholesterol"))
}
for (e in selected_exposures_bio_RBC) {
  g <- add_edges(g, c(e, "nucleated_red_blood_cell_count"))
}

# Define colors for each category
V(g)$color <- ifelse(V(g)$category == "confounder", "lightgray",
                     ifelse(grepl("exposure", V(g)$category),"lightsalmon",
                            ifelse(V(g)$category == "biomarker", "lightpink", "lightblue")))

# Set node size and label properties
V(g)$size <- 15  # default size
V(g)$size[1:(length(confounders)+length(outcome)+length(selected_bio))] <- 25


V(g)$label <- V(g)$name
V(g)$label.cex <- 0.8
V(g)$label.cex[1:(length(confounders)+length(outcome)+length(selected_bio))] <- 0.95
V(g)$label.color <- "black"
V(g)$label.font <- 4  # Options: 1 = plain, 2 = bold, 3 = italic, 4 = bold italic


# Set edge properties
E(g)$arrow.size <- 0.5
E(g)$width <- 1
E(g)$color <- "gray30"

#Change labels names
pretty_labels <- c(
  "sleep_pattern"      = "Sleep\nPattern",
  "employment_status"  = "Employment\nStatus",
  "edu"                = "Education\nLevel",
  "comp_body_10"       = "Body\nSize",
  "comp_height_10"  = "Body\nHeight",
  "social_isolation"   = "Social\nIsolation",
  "iibs_2yr"           = "IBS\n(2-Year)",
  "dist_greenspace" = "Distance\nGreen space",
  "own_rent" = "Own/Rent",
  "HDL_cholesterol" = "HDL\nCholesterol",
  "nucleated_red_blood_cell_count" = "Nucleated\nRBC Count",
  "age" = "Age",
  "ethnicity" = "Ethnicity",
  "sex"= "Sex",
  "ozone"="Ozone",
  "bmi" = "BMI",
  "imd_quin" = 	"Multiple\nDeprivation",
  "smoke_status"="Smoking\nStatus",
  "alcohol_category"="Alcohol\nCategory",
  "sys_bp" = "Systolic B.P",
  "waist_circumference"= "Waist\nCircumference",
  "met_score" ="MET\nscore",
  "physician_natd"= "Seen GP/psy\nfor anx "
)

exposures <- c(only_direct_exposures,HDL_outcome,only_HDL,all,RBC_outcome,RBD_HDL,only_RBC)

V(g)$name <- ifelse(V(g)$name %in% names(pretty_labels),
                     pretty_labels[V(g)$name],
                    V(g)$name)  # fallback to original if no pretty label

#E(g)[ .to(pretty_labels[selected_bio[2]])]$color <- "blue"
#E(g)[ .to(pretty_labels[selected_bio[1]])]$color <- "blue"

#E(g)[.from(pretty_labels[selected_exposures_outcome[!selected_exposures_outcome %in% only_direct_exposures]]) & .to(outcome)]$color <- "red"
#E(g)[.from(pretty_labels[exposures[exposures %in% RBC_outcome]]) & .to(pretty_labels[selected_bio[2]])]$color <- "grey30"
#E(g)[.from(pretty_labels[exposures[exposures %in% HDL_outcome]]) & .to(pretty_labels[selected_bio[1]])]$color <- "grey30"
#E(g)[.from(pretty_labels[confounders]) & .to(pretty_labels[selected_bio[1]])]$color <- "grey30"
#E(g)[.from(pretty_labels[confounders]) & .to(pretty_labels[selected_bio[2]])]$color <- "grey30"

#E(g)[.from(pretty_labels[only_direct_exposures])]$color <- "red"
#E(g)[.from(pretty_labels[selected_bio]) & .to(outcome) ]$color <- "blue"

# Create a custom layout that groups nodes by category

custom_grouped_layout <- function(graph) {
  n <- vcount(graph)
  layout <- matrix(0, n, 2)
  
  # Get indices by category
  confounder_indices <- which(V(graph)$category == "confounder")
  direct_exposure_indices <- which(V(graph)$category == "direct_exposure")
  indirect_exposure_indices <- which(V(graph)$category == "indirect_exposure")
  mixed_exposure_indices <- which(V(graph)$category == "mixed_exposure")
  biomarker_indices <- which(V(graph)$category == "biomarker")
  outcome_indices <- which(V(graph)$category == "outcome")
  
  # Combined exposure indices for unified y-spacing
  exposure_indices_all <- c(direct_exposure_indices, mixed_exposure_indices, indirect_exposure_indices)
  
  # Layout spacings
  x_spacing <- 2
  y_spacing <- 30
  y_spacing_exp <- 30
  
  ## Confounders - top middle
  layout[confounder_indices, 1] <- seq(-x_spacing, x_spacing, length.out = length(confounder_indices))
  layout[confounder_indices, 2] <- 500 + 3 * y_spacing
  
  ## Exposures - far left, vertically spaced
  exp_y_positions <- seq(0, y_spacing_exp * (length(exposure_indices_all) - 1), length.out = length(exposure_indices_all))
  
  layout[direct_exposure_indices, 1] <- -3 * x_spacing
  layout[direct_exposure_indices, 2] <- exp_y_positions[1:length(direct_exposure_indices)]
  
  layout[mixed_exposure_indices, 1] <- -4 * x_spacing
  layout[mixed_exposure_indices, 2] <- 10 + exp_y_positions[(length(direct_exposure_indices)+1):(length(direct_exposure_indices)+length(mixed_exposure_indices))]
  
  layout[indirect_exposure_indices, 1] <- -5 * x_spacing
  layout[indirect_exposure_indices, 2] <- 10 + exp_y_positions[(length(direct_exposure_indices)+length(mixed_exposure_indices)+1):length(exp_y_positions)]
  
  ## Biomarkers - middle
  layout[biomarker_indices, 1] <- 0
  layout[biomarker_indices, 2] <- 300 + seq(0, 5*y_spacing * (length(biomarker_indices) - 1), length.out = length(biomarker_indices))
  
  ## Outcome - right side, centered
  layout[outcome_indices, 1] <- 3 * x_spacing
  layout[outcome_indices, 2] <- mean(layout[biomarker_indices, 2])-300  # Align vertically with biomarkers
  
  return(layout)
}

# Generate layout
layout <- custom_grouped_layout(g)

#png("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/stability_results/dag_mediation_analysis.png", width = 1500, height = 1800, res = 200)

# Plot the graph
par(mar = c(1, 1, 1, 1))
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     edge.curved = 0.1,  # Slightly curved edges for better visibility
     main = "DAG for Mediation Analysis")

# Add a legend
legend("right", inset = c(-0.1, 0),
       legend = c("Confounders", "Exposures", "Biomarkers", "Outcome"),
       fill = c("lightgray", "lightsalmon", "lightpink", "lightblue"),
       cex = 1, xpd = T)

# Close the device to write the file

#dev.off()


