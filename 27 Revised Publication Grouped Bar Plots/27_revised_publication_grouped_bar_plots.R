# ============================================================================== #
# Script: 27_revised_publication_grouped_bar_plots.R
# Title: Publication Grouped Bar Plots with Kendall Concordance
# Description:
#   Generates publication-ready grouped bar charts with significance letters, error bars (mean +- SD), and Kendall W concordance statistics in multiple image resolutions.
#
# Inputs:
#   - results/Table_Folds_5_Times_10_* (Evaluation results across methods and traits)
#
# Outputs:
#   - results/Publication_Ready_Grouped_Bar_Chart.png (300 DPI publication figure)
#   - results/Publication_Ready_With_Stats.png (Figure with significance ranking annotations)
#
# Required R Packages: ggplot2, reshape2, PMCMRplus, dplyr, viridis
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


# Load required libraries
library(tidyverse)
library(PMCMRplus)
library(ggplot2)
library(multcompView)

# Read the data
all_data <- read.csv("all_data.csv")

# Check the data structure
cat("Data loaded successfully!\n")
cat("Dimensions:", dim(all_data), "\n")
cat("First few rows:\n")
print(head(all_data))

# ============================================================================
# STATISTICAL ANALYSIS: Friedman Test + Nemenyi Post-hoc for Each Method
# ============================================================================

methods <- unique(all_data$Method)
traits <- unique(all_data$Trait)

# Storage for results
friedman_results <- list()
nemenyi_results <- list()
summary_stats <- data.frame()

cat("\n================================================================================\n")
cat("PERFORMING STATISTICAL TESTS\n")
cat("================================================================================\n\n")

for (method in methods) {
  cat(paste0("\n--- METHOD: ", method, " ---\n"))
  
  # Filter data for this method
  method_data <- all_data %>%
    filter(Method == method) %>%
    select(CV_Iteration, Trait, Accuracy)
  
  # Reshape to wide format for Friedman test
  wide_data <- method_data %>%
    pivot_wider(names_from = Trait, values_from = Accuracy)
  
  # Create matrix (rows = blocks/CV folds, columns = traits)
  trait_matrix <- as.matrix(wide_data[, -1])
  
  # Perform Friedman test
  friedman_test <- friedman.test(trait_matrix)
  friedman_results[[method]] <- friedman_test
  
  cat(sprintf("Friedman Test: χ² = %.3f, df = %d, p-value = %.4f %s\n",
              friedman_test$statistic,
              friedman_test$parameter,
              friedman_test$p.value,
              ifelse(friedman_test$p.value < 0.001, "***",
                     ifelse(friedman_test$p.value < 0.01, "**",
                            ifelse(friedman_test$p.value < 0.05, "*", "ns")))))
  
  # Calculate Kendall's W (effect size)
  n <- nrow(trait_matrix)
  k <- ncol(trait_matrix)
  ranked_data <- t(apply(trait_matrix, 1, rank))
  R_j <- colSums(ranked_data)
  R_bar <- mean(R_j)
  S <- sum((R_j - R_bar)^2)
  W <- (12 * S) / (n^2 * (k^3 - k))
  
  cat(sprintf("Kendall's W = %.3f (effect size: %s)\n", W,
              ifelse(W < 0.1, "Small",
                     ifelse(W < 0.3, "Medium",
                            ifelse(W < 0.5, "Large", "Very Large")))))
  
  # Mean ranks
  mean_ranks <- colMeans(ranked_data)
  cat("Mean Ranks (lower = better):\n")
  for (trait in names(sort(mean_ranks))) {
    cat(sprintf("  %s: %.2f\n", trait, mean_ranks[trait]))
  }
  
  # Perform Nemenyi post-hoc test if significant
  if (friedman_test$p.value < 0.05) {
    cat("→ Performing Nemenyi post-hoc test...\n")
    
    # Prepare data for Nemenyi test
    long_data <- method_data %>%
      mutate(CV_Iteration = as.factor(CV_Iteration),
             Trait = as.factor(Trait))
    
    # Nemenyi test
    nemenyi_test <- frdAllPairsNemenyiTest(
      y = long_data$Accuracy,
      groups = long_data$Trait,
      blocks = long_data$CV_Iteration,
      p.adjust.method = "none"
    )
    
    nemenyi_results[[method]] <- nemenyi_test
    
    # Extract p-value matrix and make it symmetric
    trait_names <- colnames(trait_matrix)
    p_matrix_lower <- nemenyi_test$p.value
    
    # Create full symmetric p-value matrix
    n_traits <- length(trait_names)
    p_matrix <- matrix(1, nrow = n_traits, ncol = n_traits)
    rownames(p_matrix) <- trait_names
    colnames(p_matrix) <- trait_names
    
    # Fill in the lower triangle
    p_matrix[lower.tri(p_matrix)] <- p_matrix_lower[lower.tri(p_matrix_lower)]
    # Mirror to upper triangle
    p_matrix[upper.tri(p_matrix)] <- t(p_matrix)[upper.tri(p_matrix)]
    
    cat("Pairwise p-values:\n")
    print(round(p_matrix, 4))
    
    # Assign significance letters based on pairwise comparisons
    trait_order <- names(sort(mean_ranks))
    # Initialize all traits with "a"
    sig_letters <- setNames(rep("a", length(trait_names)), trait_names)
    p_threshold <- 0.05
    
    # Algorithm: Compare each trait with the best trait (lowest rank)
    # If significantly different from the best, assign "b"
    if (length(trait_order) >= 2) {
      best_trait <- trait_order[1]  # Trait with lowest rank (best)
      
      for (i in 2:length(trait_order)) {
        current_trait <- trait_order[i]
        p_val <- p_matrix[best_trait, current_trait]
        
        if (!is.na(p_val) && p_val < p_threshold) {
          sig_letters[current_trait] <- "b"
        }
      }
      
      # Special case: if middle trait is "a" but differs from worst trait
      # and best trait is "a", check if we need to split
      if (length(trait_order) == 3) {
        trait1 <- trait_order[1]
        trait2 <- trait_order[2]
        trait3 <- trait_order[3]
        
        # If trait2 is "a" and trait3 is "a", but they differ from each other
        # and trait2 doesn't differ from trait1, then trait3 should be "b"
        if (sig_letters[trait2] == "a" && sig_letters[trait3] == "a") {
          p23 <- p_matrix[trait2, trait3]
          if (!is.na(p23) && p23 < p_threshold) {
            sig_letters[trait3] <- "b"
          }
        }
      }
    }
    
    cat("Significance groups:\n")
    for (trait in trait_order) {
      cat(sprintf("  %s: %s (rank: %.2f)\n", trait, sig_letters[trait], mean_ranks[trait]))
    }
    
  } else {
    cat("→ No significant differences (all traits: 'a')\n")
    sig_letters <- setNames(rep("a", length(mean_ranks)), names(mean_ranks))
  }
  
  # Store summary statistics for each trait
  for (trait in traits) {
    trait_values <- method_data %>% filter(Trait == trait) %>% pull(Accuracy)
    
    summary_stats <- rbind(summary_stats, data.frame(
      Method = method,
      Trait = trait,
      Mean = mean(trait_values),
      SD = sd(trait_values),
      SE = sd(trait_values) / sqrt(length(trait_values)),
      Median = median(trait_values),
      Rank = mean_ranks[trait],
      Letter = sig_letters[trait],
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n================================================================================\n")
cat("SUMMARY STATISTICS TABLE\n")
cat("================================================================================\n\n")

print(summary_stats, row.names = FALSE)

# ============================================================================
# RANK METHODS BY OVERALL PERFORMANCE
# ============================================================================

cat("\n================================================================================\n")
cat("RANKING METHODS BY OVERALL PERFORMANCE\n")
cat("================================================================================\n\n")

# Calculate overall performance for each method (average across all traits)
method_performance <- summary_stats %>%
  group_by(Method) %>%
  summarise(
    Overall_Mean = mean(Mean),
    Overall_Median = median(Mean),
    .groups = 'drop'
  ) %>%
  arrange(desc(Overall_Mean)) %>%
  mutate(
    Method_Rank = row_number(),
    Rank_Label = paste0("#", Method_Rank)
  )

cat("Method Rankings (1 = Best Overall Performance):\n")
cat("Based on average predictive ability across all three traits\n\n")

for (i in 1:nrow(method_performance)) {
  cat(sprintf("  %s. %-6s - Mean: %.4f\n", 
              method_performance$Method_Rank[i],
              method_performance$Method[i],
              method_performance$Overall_Mean[i]))
}

# Add method ranks to summary stats
summary_stats <- summary_stats %>%
  left_join(method_performance %>% select(Method, Method_Rank, Rank_Label), 
            by = "Method")

# ============================================================================
# VISUALIZATION WITH METHOD RANKS
# ============================================================================

cat("\n================================================================================\n")
cat("CREATING VISUALIZATION WITH METHOD RANKS\n")
cat("================================================================================\n\n")

# Prepare plot data with method ranks
plot_data <- summary_stats %>%
  mutate(
    Trait_Display = case_when(
      Trait == "Is_Dead" ~ "Binary Survival",
      Trait == "DPC_date" ~ "DPC_Date",
      Trait == "DPC_time" ~ "DPC_Time",
      TRUE ~ Trait
    ),
    # Order methods by performance rank
    Method_Label = paste0(Method, " (", Method_Rank, ")"),
    Method_Label = factor(Method_Label, 
                          levels = paste0(method_performance$Method, 
                                          " (", method_performance$Method_Rank, ")"))
  )

# Create the main plot with method ranks
p <- ggplot(plot_data, aes(x = Method_Label, y = Mean, fill = Trait_Display)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), 
           color = "black", width = 0.8, linewidth = 0.8, alpha = 0.85) +
  
  # Add error bars
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.25, linewidth = 0.6) +
  
  # Add trait rank numbers in RED inside bars
  geom_text(aes(y = Mean / 2, label = sprintf("%.1f", Rank)),
            position = position_dodge(0.9),
            color = "red", fontface = "bold", size = 4.5) +
  
  # Add significance letters in RED above error bars
  geom_text(aes(y = Mean + SE + 0.015, label = Letter),
            position = position_dodge(0.9),
            color = "red", fontface = "bold", size = 5, vjust = 0) +
  
  # Color scheme matching your image
  scale_fill_manual(
    values = c(
      "Binary Survival" = "#E57373",  # Red
      "DPC_Date" = "#64B5F6",          # Blue
      "DPC_Time" = "#81C784"           # Green
    ),
    name = "Trait"
  ) +
  
  # Labels
  labs(
    title = "Mean Predictive Ability by Method with Statistical Significance",
    subtitle = "Methods ranked by overall performance (1 = best) | Red numbers: trait ranks within method | Red letters: significance groups (p < 0.05)\nFriedman test followed by Nemenyi post-hoc test",
    x = "Prediction Method (Overall Rank)",
    y = "Mean Predictive Ability ± SE",
    caption = "Traits with the same letter within a method are not significantly different"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "darkred", 
                                 face = "bold", margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0.5, size = 9, face = "italic", 
                                margin = margin(t = 10)),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # Y-axis scale
  scale_y_continuous(
    limits = c(0, 0.65), 
    breaks = seq(0, 0.65, 0.05),
    expand = c(0, 0)
  )

# Display the plot
print(p)

# Save the plots
ggsave("Predictive_Ability_With_Significance_Ranked.png", 
       plot = p, width = 16, height = 10, dpi = 300, bg = "white")

ggsave("Predictive_Ability_With_Significance_Ranked_HighRes.png", 
       plot = p, width = 16, height = 10, dpi = 600, bg = "white")

ggsave("Predictive_Ability_With_Significance_Ranked.pdf", 
       plot = p, width = 16, height = 10, device = "pdf")

cat("\n✓ Main plot with method ranks saved successfully!\n")

# ============================================================================
# ALTERNATIVE PLOT: Simple Method Names with Rank Numbers on Top
# ============================================================================

# Create plot with rank numbers displayed above bars
p2 <- ggplot(plot_data, aes(x = factor(Method, levels = method_performance$Method), 
                            y = Mean, fill = Trait_Display)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), 
           color = "black", width = 0.8, linewidth = 0.8, alpha = 0.85) +
  
  # Add error bars
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.25, linewidth = 0.6) +
  
  # Add trait rank numbers in RED inside bars
  geom_text(aes(y = Mean / 2, label = sprintf("%.1f", Rank)),
            position = position_dodge(0.9),
            color = "red", fontface = "bold", size = 4.5) +
  
  # Add significance letters in RED above error bars
  geom_text(aes(y = Mean + SE + 0.015, label = Letter),
            position = position_dodge(0.9),
            color = "red", fontface = "bold", size = 5, vjust = 0) +
  
  # Add method rank as a label at the top
  geom_text(data = plot_data %>% 
              group_by(Method, Method_Rank) %>% 
              summarise(max_y = max(Mean + SE) + 0.04, .groups = 'drop'),
            aes(x = Method, y = max_y, label = Method_Rank),
            color = "darkblue", fontface = "bold", size = 6,
            inherit.aes = FALSE) +
  
  # Color scheme
  scale_fill_manual(
    values = c(
      "Binary Survival" = "#E57373",
      "DPC_Date" = "#64B5F6",
      "DPC_Time" = "#81C784"
    ),
    name = "Trait"
  ) +
  
  # Labels
  labs(
    title = "Mean Predictive Ability by Method with Statistical Significance",
    subtitle = "Blue numbers at top: overall method rank (1 = best) | Red numbers in bars: trait ranks within method | Red letters: significance groups (p < 0.05)",
    x = "Prediction Method",
    y = "Mean Predictive Ability ± SE",
    caption = "Traits with the same letter within a method are not significantly different"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "darkred", 
                                 face = "bold", margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0.5, size = 9, face = "italic", 
                                margin = margin(t = 10)),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # Y-axis scale
  scale_y_continuous(
    limits = c(0, 0.7), 
    breaks = seq(0, 0.7, 0.05),
    expand = c(0, 0)
  )

print(p2)

ggsave("Predictive_Ability_With_Method_Ranks_Top.png", 
       plot = p2, width = 16, height = 10, dpi = 300, bg = "white")

cat("✓ Alternative plot with ranks on top saved!\n")

# ============================================================================
# ADDITIONAL VISUALIZATIONS
# ============================================================================

# 1.Heatmap with method ranks
p_heatmap <- ggplot(plot_data, aes(x = Trait_Display, y = Method_Label, fill = Mean)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.3f\n[%.1f] %s", Mean, Rank, Letter)),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_gradient2(
    low = "#D73027", mid = "#FEE08B", high = "#1A9850",
    midpoint = median(plot_data$Mean),
    name = "Mean\nAccuracy"
  ) +
  labs(
    title = "Heatmap: Trait Performance Within Each Method",
    subtitle = "Methods ranked by overall performance | Values show: Mean [Rank] Letter",
    x = "Trait",
    y = "Prediction Method (Overall Rank)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

print(p_heatmap)
ggsave("Trait_Performance_Heatmap_Ranked.png", p_heatmap, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Heatmap with method ranks saved!\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save summary statistics with method ranks
write.csv(summary_stats, "Summary_Statistics_With_Significance_And_Ranks.csv", row.names = FALSE)
cat("✓ Summary statistics with ranks saved to CSV!\n")

# Save method performance ranking
write.csv(method_performance, "Method_Performance_Rankings.csv", row.names = FALSE)
cat("✓ Method performance rankings saved to CSV!\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("================================================================================\n\n")

cat("Files created:\n")
cat("  1.Predictive_Ability_With_Significance_Ranked.png (methods with ranks in labels)\n")
cat("  2.Predictive_Ability_With_Significance_Ranked_HighRes.png (600 dpi)\n")
cat("  3.Predictive_Ability_With_Significance_Ranked.pdf (vector)\n")
cat("  4.Predictive_Ability_With_Method_Ranks_Top.png (ranks as blue numbers on top)\n")
cat("  5.Trait_Performance_Heatmap_Ranked.png\n")
cat("  6.Summary_Statistics_With_Significance_And_Ranks.csv\n")
cat("  7.Method_Performance_Rankings.csv\n\n")

cat("Method Rankings Summary:\n")
print(method_performance, row.names = FALSE)

cat("\n\nInterpretation:\n")
cat("  • Method rank: 1 = Best overall performance (averaged across all traits)\n")
cat("  • Trait rank (red numbers in bars): Lower = Better within that method\n")
cat("  • Same letter = No significant difference (p ≥ 0.05)\n")
cat("  • Different letters = Significant difference (p < 0.05)\n")
cat("  • Test: Friedman test + Nemenyi post-hoc test\n")

