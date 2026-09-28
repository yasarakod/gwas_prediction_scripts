# ============================================================================== #
# Script: 18_gblup_popsize_folds_trajectories.R
# Title: GBLUP Population Size & Fold Trajectory Plots
# Description:
#   Generates trajectory line plots showing mean predictive ability and accuracy trends across population sizes for 3-fold vs 5-fold cross-validation.
#
# Inputs:
#   - results/GBLUP_Folds_3_Times_10_* and GBLUP_Folds_5_Times_10_* (Output result files)
#
# Outputs:
#   - results/Table_GBLUP_Combined.csv (Combined accuracy dataframe)
#   - results/Plot_GBLUP_Combined.png (Multi-line trajectory plot for accuracy)
#   - results/Ability_Plot_GBLUP_Combined.png (Multi-line trajectory plot for ability)
#
# Required R Packages: ggplot2, reshape2, dplyr
# ============================================================================== #

# ------------------------------------------------------------------------------
# Directory Setup & Output Path Configuration
# ------------------------------------------------------------------------------
# Automatically create "results" folder for output deliverables
results_dir <- file.path(getwd(), "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

# Automatically locate shared input data directory
if (dir.exists(file.path(getwd(), "data"))) {
  data_dir <- file.path(getwd(), "data")
} else if (dir.exists(file.path(dirname(getwd()), "data"))) {
  data_dir <- file.path(dirname(getwd()), "data")
} else {
  data_dir <- getwd()
}
data_path <- paste0(data_dir, "/")
# ------------------------------------------------------------------------------


############################After running CV code###############################
###series of plots for different population sizes and cross validation folds####
# setwd() configured via project environment # Set working directory

# Parse Arguments
dir<-"data"
graph_type<-"C"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("DPC_time")
pmethods<-c("GBLUP")
iter<-10
fold<-c(3,5)
#fold<-5
pop_sizes = c("50","100", "150", "200", "250", "300", "350", "400", "474")
marker_sizes = c("top_1000")
#marker_sizes<-c("top_40000")
#marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("NULL")
h2<-1

library(reshape2)
library(ggplot2)

data_path <- paste0(dir, "/")  # make sure folder path ends with /

if (graph_type == "C") {
  
  sname <- trait_name
  pmethod <- pmethods
  
  # Initialize an empty data.frame to store all results
  all_data <- data.frame()
  
  for (fold_val in fold) {
    for (i in seq_along(pop_sizes)) {
      
      # Build filename
      filename <- paste0(
        data_path, pmethod, "_Folds_", fold_val, "_Times_", iter, "_",
        sname, "_", Fixed_effect, "_", pop_sizes[i], "_", marker_sizes, ".txt"
      )
      
      # Read predictions
      datap <- read.table(filename)
      
      # Compute correlation for each iteration
      for (n in 1:iter) {
        temp <- data.frame(
          Fold = fold_val,
          PopSize = pop_sizes[i],
          Iteration = n,
          Correlation = cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
        )
        all_data <- rbind(all_data, temp)
      }
    }
  }
  
  # Optional: normalize by heritability
  all_data$Correlation_norm <- all_data$Correlation / sqrt(h2)
  
  # Save combined CSV
  write.csv(all_data, paste0(data_path, "Table_", pmethod, "_Combined.csv"), row.names = FALSE)
  
  # Plot: normalized correlation
  p <- ggplot(all_data, aes(x = factor(PopSize), y = Correlation_norm, fill = factor(PopSize))) +
    geom_boxplot(color = "black") +
    labs(title = paste0("Prediction Accuracy: ", sname),
         x = "Population Size", y = "Prediction Accuracy (normalized)",
         fill = "Pop Size") +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "red") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(paste0(data_path, "Plot_", pmethod, "_Combined.png"), p, dpi = 300, width = 10, height = 6)
  
  # Plot: raw correlation (without heritability normalization)
  p2 <- ggplot(all_data, aes(x = factor(PopSize), y = Correlation, fill = factor(PopSize))) +
    geom_boxplot(color = "black") +
    labs(title = paste0("Predictive Ability: ", sname),
         x = "Population Size", y = "Predictive Ability",
         fill = "Pop Size") +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "red") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(paste0(data_path, "Ability_Plot_", pmethod, "_Combined.png"), p2, dpi = 300, width = 10, height = 6)
}
#save all_data
write.csv(all_data, paste0(data_path, "Table_", pmethod, "_Combined.csv"), row.names = FALSE)

library(ggplot2)

# Define desired order: 50 first, then ascending
desired_order <- c(50, 100, 150, 200, 250, 300, 350, 400, 474)

# Convert PopSize to a factor with the specified order
all_data$PopSize <- factor(all_data$PopSize, levels = desired_order, ordered = TRUE)

# Fix the Fold levels to match legend labels
all_data$Fold <- factor(all_data$Fold, levels = c("3", "5"), labels = c("3 x 10", "5 x 10"))

# Create boxplot with separate colors for each Fold
p <- ggplot(all_data, aes(x = PopSize, y = Correlation, fill = Fold)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  stat_summary(aes(group = Fold),
               fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  labs(title = "GBLUP_DPC_time",
       x = "Population Size", 
       y = "Predictive Ability",
       fill = "Cross Validation") +
  scale_fill_manual(values = c("3 x 10" = "#1F78B4", "5 x 10" = "#33A02C")) +
  theme_minimal() +
  theme(
    legend.title = element_text(colour = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  )

print(p)
ggsave(filename = "GBLUP_PopSize_Folds_DPC_date.png", plot = p, width = 10, height = 6, dpi = 300)

library(ggplot2)

# Define desired order: 50 first, then ascending
desired_order <- c(50, 100, 150, 200, 250, 300, 350, 400, 474)
all_data$PopSize <- factor(all_data$PopSize, levels = desired_order, ordered = TRUE)
all_data$Fold <- factor(all_data$Fold, levels = c("3", "5"), labels = c("3 x 10", "5 x 10"))

# Calculate mean Correlation for each PopSize and Fold
mean_data <- aggregate(Correlation ~ PopSize + Fold, all_data, mean)

p <- ggplot(all_data, aes(x = PopSize, y = Correlation, fill = Fold)) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  stat_summary(aes(group = Fold),
               fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  # Add dashed lines connecting the means for each Fold
  geom_line(
    data = mean_data,
    aes(x = PopSize, y = Correlation, group = Fold, color = Fold),
    linetype = "dashed", linewidth = 1.2, show.legend = FALSE
  ) +
  scale_color_manual(values = c("3 x 10" = "#1F78B4", "5 x 10" = "#33A02C")) +
  labs(title = "GBLUP_DPC_time",
       x = "Population Size", 
       y = "Predictive Ability",
       fill = "Cross Validation") +
  scale_fill_manual(values = c("3 x 10" = "#1F78B4", "5 x 10" = "#33A02C")) +
  theme_minimal() +
  theme(
    legend.title = element_text(colour = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  )
print(p)
ggsave(filename = "GBLUP_PopSize_Folds_DPC_time_MeanLines.png", plot = p, width = 10, height = 6, dpi = 300)
