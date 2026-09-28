# ============================================================================== #
# Script: 20_friedman_nemenyi_analysis.R
# Title: Friedman & Nemenyi Non-Parametric Model Ranking
# Description:
#   Executes Friedman two-way ANOVA by ranks and Nemenyi post-hoc tests on cross-validation accuracy matrices, outputting statistical rankings and p-value tables.
#
# Inputs:
#   - results/Table_Folds_5_Times_10_*_top_40000.csv (Consolidated accuracy tables across prediction methods)
#
# Outputs:
#   - results/Friedman_Results_*.csv (Friedman test statistics)
#   - results/Method_Rankings_*.csv (Method rankings and significance groups)
#   - results/Plot_with_Rankings_*.png (Boxplot annotated with method rankings)
#
# Required R Packages: tidyverse, PMCMRplus, rstatix, ggplot2
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


##Once the cross-validation is done, After the results are stored in a text file#####
##########The following code can be used to generate the graphs and tables#####
# ---- Argument Handling (Combined Approach) ----
# setwd() configured via project environment

# Load required libraries
library(tidyverse)
library(PMCMRplus)
install.packages("PMCMRplus", dependencies = TRUE)
library(rstatix)
library(ggplot2)
library(reshape2)

# Read the CSV file
data <- read.csv("Table_Folds_5_Times_10_DPC_date_NULL_NULL_top_40000.csv", row.names = 1)

# Display the data structure
cat("Data dimensions:", dim(data), "\n")
cat("Prediction methods:", colnames(data), "\n")
cat("Cross-validation folds:", rownames(data), "\n\n")

# Convert data to long format for analysis
data_long <- data %>%
  rownames_to_column("CV_Fold") %>%
  pivot_longer(cols = -CV_Fold, 
               names_to = "Method", 
               values_to = "Prediction_Ability")

# Display summary statistics
cat("=== SUMMARY STATISTICS ===\n")
summary_stats <- data_long %>%
  group_by(Method) %>%
  summarise(
    Mean = round(mean(Prediction_Ability), 4),
    SD = round(sd(Prediction_Ability), 4),
    Median = round(median(Prediction_Ability), 4),
    Min = round(min(Prediction_Ability), 4),
    Max = round(max(Prediction_Ability), 4),
    .groups = 'drop'
  ) %>%
  arrange(desc(Mean))

print(summary_stats)

# Create a boxplot for visual comparison
cat("\n=== CREATING BOXPLOT ===\n")
p1 <- ggplot(data_long, aes(x = reorder(Method, Prediction_Ability, median), 
                           y = Prediction_Ability, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Prediction Ability by Method",
       x = "Prediction Method",
       y = "Prediction Ability") +
  coord_flip()

print(p1)

# Perform Friedman test
cat("\n=== FRIEDMAN TEST ===\n")
# Prepare data matrix for Friedman test (rows = blocks/CV folds, columns = treatments/methods)
data_matrix <- as.matrix(data)

# Friedman test
friedman_result <- friedman.test(data_matrix)
print(friedman_result)

# Calculate effect size (Kendall's W)
n_blocks <- nrow(data_matrix)
n_treatments <- ncol(data_matrix)
chi_squared <- friedman_result$statistic
kendall_w <- chi_squared / (n_blocks * (n_treatments - 1))

cat("\nEffect size (Kendall's W):", round(kendall_w, 4))
cat("\nInterpretation of Kendall's W:")
cat("\n  - Small effect: W < 0.1")
cat("\n  - Medium effect: 0.1 ≤ W < 0.3") 
cat("\n  - Large effect: W ≥ 0.3")

if(kendall_w < 0.1) {
  cat("\n  -> Small effect size")
} else if(kendall_w < 0.3) {
  cat("\n  -> Medium effect size")
} else {
  cat("\n  -> Large effect size")
}

# Post-hoc analysis if Friedman test is significant
if(friedman_result$p.value < 0.05) {
  cat("\n\n=== POST-HOC ANALYSIS ===")
  cat("\nFriedman test is significant (p < 0.05). Proceeding with post-hoc tests...\n")
  
  # Nemenyi post-hoc test (recommended for Friedman test)
  cat("\n--- Nemenyi Post-hoc Test ---\n")
  nemenyi_result <- PMCMRplus::frdAllPairsNemenyiTest(data_matrix)
  print(nemenyi_result)
  
  # Create a summary of significant pairwise comparisons
  cat("\n--- Summary of Significant Pairwise Comparisons (p < 0.05) ---\n")
  p_values <- nemenyi_result$p.value
  significant_pairs <- which(p_values < 0.05, arr.ind = TRUE)
  
  if(nrow(significant_pairs) > 0) {
    for(i in 1:nrow(significant_pairs)) {
      row_idx <- significant_pairs[i, 1]
      col_idx <- significant_pairs[i, 2]
      method1 <- rownames(p_values)[row_idx]
      method2 <- colnames(p_values)[col_idx]
      p_val <- p_values[row_idx, col_idx]
      cat(sprintf("%s vs %s: p = %.4f\n", method1, method2, p_val))
    }
  } else {
    cat("No significant pairwise differences found in post-hoc analysis.\n")
  }
  
  # Alternative: Conover post-hoc test (more powerful but less conservative)
  cat("\n--- Conover Post-hoc Test (Alternative) ---\n")
  conover_result <- PMCMRplus::frdAllPairsConoverTest(data_matrix)
  print(conover_result)
  
  # Rank-based analysis
  cat("\n--- Mean Ranks ---\n")
  ranks <- apply(data_matrix, 1, rank)
  mean_ranks <- rowMeans(ranks)
  rank_summary <- data.frame(
    Method = names(mean_ranks),
    Mean_Rank = round(mean_ranks, 2)
  ) %>%
  arrange(Mean_Rank)
  
  print(rank_summary)
  
} else {
  cat("\n\nFriedman test is not significant (p ≥ 0.05).")
  cat("\nNo significant differences found among prediction methods.")
  cat("\nPost-hoc analysis is not recommended.\n")
}

# Additional visualization: Mean ranks plot
cat("\n=== CREATING MEAN RANKS PLOT ===\n")
if(exists("mean_ranks")) {
  rank_df <- data.frame(
    Method = names(mean_ranks),
    Mean_Rank = mean_ranks
  )
  
  p2 <- ggplot(rank_df, aes(x = reorder(Method, -Mean_Rank), y = Mean_Rank, fill = Method)) +
    geom_col(alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Mean Ranks of Prediction Methods",
         subtitle = "Lower ranks indicate better performance",
         x = "Prediction Method",
         y = "Mean Rank") +
    geom_text(aes(label = round(Mean_Rank, 2)), vjust = -0.5)
  
  print(p2)
}

# Save results
cat("\n=== SAVING RESULTS ===\n")
# Save summary statistics
write.csv(summary_stats, "prediction_methods_summary.csv", row.names = FALSE)

# Save detailed results
results_summary <- list(
  friedman_test = friedman_result,
  effect_size_kendall_w = kendall_w,
  summary_statistics = summary_stats
)

if(exists("nemenyi_result")) {
  results_summary$nemenyi_posthoc <- nemenyi_result
  results_summary$mean_ranks <- rank_summary
}

# Save plots
ggsave("prediction_methods_boxplot.png", p1, width = 10, height = 8, dpi = 300)
if(exists("p2")) {
  ggsave("mean_ranks_plot.png", p2, width = 10, height = 6, dpi = 300)
}

cat("Analysis complete!")
cat("\nFiles saved:")
cat("\n- prediction_methods_summary.csv")
cat("\n- prediction_methods_boxplot.png")
if(exists("p2")) {
  cat("\n- mean_ranks_plot.png")
}

# Print interpretation guidelines
cat("\n\n=== INTERPRETATION GUIDELINES ===")
cat("\n1.Friedman Test:")
cat("\n   - Tests if there are significant differences among methods")
cat("\n   - H0: All methods perform equally")
cat("\n   - H1: At least one method differs significantly")

cat("\n\n2.Effect Size (Kendall's W):")
cat("\n   - Measures the strength of agreement among CV folds")
cat("\n   - Range: 0 (no agreement) to 1 (perfect agreement)")

cat("\n\n3.Post-hoc Tests (if Friedman is significant):")
cat("\n   - Nemenyi: Conservative, controls family-wise error rate")
cat("\n   - Conover: More powerful, may detect more differences")

cat("\n\n4.Mean Ranks:")
cat("\n   - Lower rank = better performance")
cat("\n   - Methods with significantly different ranks perform differently")

cat("\n\nNote: In your data, higher prediction ability values appear to be better,")
cat("\nso methods with lower mean ranks are actually performing better.\n")
#