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

################### FULL GRAPH ########################### 


#load data 

X_biomarkers <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_internal.rds"))
X_exposures <- as.data.frame(readRDS("/rds/general/project/hda_24-25/live/Comp_Epi/Group03/UKB_data/data_external.rds"))
biomarkers <- colnames(X_biomarkers)
exposures <- colnames(X_exposures)

confounders <- c("age", "sex", "ethnicity") 
outcome <- c("CKD")
exposures <- exposures[ !exposures %in% confounders]

# Create a list of all nodes
all_nodes <- c(confounders, exposures, biomarkers, outcome)

# Create node category indicator (for coloring and grouping)
node_category <- c(
  rep("confounder", length(confounders)),
  rep("exposure", length(exposures)),
  rep("biomarker", length(biomarkers)),
  rep("outcome", length(outcome))
)

# Define edges for the DAG - only between different categories
edges <- c()

# Create the graph with vertices but no edges
g <- make_empty_graph(n = length(all_nodes), directed = TRUE)
# Set node attributes
V(g)$name <- all_nodes
V(g)$category <- node_category

for (b in biomarkers) {
  g <- add_edges(g, c(b,outcome))
}
for (e in exposures) {
  g <- add_edges(g, c(e, outcome))
}
for (c in confounders) {
  g <- add_edges(g, c(c, outcome))
}
for (e in exposures)  {
  for (b in biomarkers){
    g <- add_edges(g, c(e, b))
  }
}

# Define colors for each category
V(g)$color <- ifelse(V(g)$category == "confounder", "lightgray",
                     ifelse(V(g)$category == "exposure", "lightsalmon",
                            ifelse(V(g)$category == "biomarker", "lightpink", "lightblue")))

# Set node size and label properties
V(g)$size <- 5
V(g)$label <- V(g)$name
V(g)$label.cex <- 0.8

# Set edge properties
E(g)$arrow.size <- 0.5
E(g)$width <- 1
E(g)$color <- "gray30"



# Create a custom layout that groups nodes by category
custom_grouped_layout <- function(graph) {
  # Get the number of nodes
  n <- vcount(graph)
  
  # Initialize layout with random positions
  layout <- matrix(0, n, 2)
  
  # Define positions for each category
  confounder_indices <- which(V(graph)$category == "confounder")
  exposure_indices <- which(V(graph)$category == "exposure")
  biomarker_indices <- which(V(graph)$category == "biomarker")
  outcome_indices <- which(V(graph)$category == "outcome")
  
  # Place confounders at the top
  confounder_x <- seq(-2, 2, length.out = length(confounder_indices))
  layout[confounder_indices, 1] <- confounder_x
  layout[confounder_indices, 2] <- 50*length(exposure_indices) + 100*length(biomarker_indices) +300
  
  # Place exposures on the left
  exposure_y <- seq(0, 100*length(exposure_indices), length.out = length(exposure_indices))
  layout[exposure_indices, 1] <- -10
  layout[exposure_indices, 2] <- exposure_y
  
  # Place biomarkers in the middle
  biomarker_y <- seq(50*length(exposure_indices),  50*length(exposure_indices) + 100*length(biomarker_indices), length.out = length(biomarker_indices))
  layout[biomarker_indices, 1] <- 0
  layout[biomarker_indices, 2] <- biomarker_y
  
  # Place outcome on the right
  layout[outcome_indices, 1] <- 10
  layout[outcome_indices, 2] <- 5
  
  return(layout)
}



# Generate layout
layout <- custom_grouped_layout(g)

# Plot the graph
par(mar = c(1, 1, 1, 1))
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     edge.curved = 0.2,  # Slightly curved edges for better visibility
     main = "DAG for Mention Analysis")

# Add a legend
legend("bottomright", 
       legend = c("Confounders", "Exposures", "Biomarkers", "outcome"),
       fill = c("lightgray", "lightsalmon", "lightpink", "lightblue"),
       cex = 0.8)


################### SPARSE GRAPH ########################### 

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

exposures <- union(union(selected_exposures_bio_HDL,selected_exposures_bio_RBC),selected_exposures_outcome)

# Create a list of all nodes
all_nodes <- c(confounders, exposures, selected_bio, outcome)
node_category <- c(
  rep("confounder", length(confounders)),
  rep("exposure", length(exposures)),
  rep("biomarker", length(selected_bio)),
  rep("outcome", length(outcome))
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
  for (e in exposures){
    g <- add_edges(g, c(c, e))
  }
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
                     ifelse(V(g)$category == "exposure", "lightsalmon",
                            ifelse(V(g)$category == "biomarker", "lightpink", "lightblue")))

# Set node size and label properties
V(g)$size <- 15
V(g)$label <- V(g)$name
V(g)$label.cex <- 0.8

# Set edge properties
E(g)$arrow.size <- 0.5
E(g)$width <- 1
E(g)$color <- "gray30"

# Create a custom layout that groups nodes by category

custom_grouped_layout <- function(graph) {
  n <- vcount(graph)
  layout <- matrix(0, n, 2)
  
  # Get indices by category
  confounder_indices <- which(V(graph)$category == "confounder")
  exposure_indices <- which(V(graph)$category == "exposure")
  biomarker_indices <- which(V(graph)$category == "biomarker")
  outcome_indices <- which(V(graph)$category == "outcome")
  
  # Layout dimensions
  x_spacing <- 2
  y_spacing <- 10
  
  # Confounders - top middle
  layout[confounder_indices, 1] <- seq(-x_spacing, x_spacing, length.out = length(confounder_indices))
  layout[confounder_indices, 2] <- 80 + 3 * y_spacing
  
  # Exposures - bottom left
  layout[exposure_indices, 1] <- -3 * x_spacing
  layout[exposure_indices, 2] <- seq(0, y_spacing * (length(exposure_indices) - 1), length.out = length(exposure_indices))
  
  # Biomarkers - middle
  layout[biomarker_indices, 1] <- 0
  layout[biomarker_indices, 2] <- 50 + seq(0, y_spacing * (length(biomarker_indices) - 1), length.out = length(biomarker_indices))
  
  # Outcome - bottom right
  layout[outcome_indices, 1] <- 3 * x_spacing
  layout[outcome_indices, 2] <- y_spacing  # Place near center vertically
  
  return(layout)
}

# Generate layout
layout <- custom_grouped_layout(g)

# Plot the graph
par(mar = c(1, 1, 1, 1))
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     edge.curved = 0.2,  # Slightly curved edges for better visibility
     main = "DAG for Mediation Analysis")

# Add a legend
legend("bottomright", 
       legend = c("Confounders", "Exposures", "Biomarkers", "outcome"),
       fill = c("lightgray", "lightsalmon", "lightpink", "lightblue"),
       cex = 0.8)
