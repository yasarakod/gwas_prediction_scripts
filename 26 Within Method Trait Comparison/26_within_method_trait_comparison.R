# ============================================================================== #
# Script: 26_within_method_trait_comparison.R
# Title: Within-Method Trait Predictability Comparison
# Description:
#   Performs statistical comparison of prediction performance across traits within each prediction model, generating performance heatmaps and trait ranking letter tables.
#
# Inputs:
#   - results/*_Folds_5_Times_10_* (Model prediction iteration results across traits)
#
# Outputs:
#   - results/Predictive_Ability_With_Significance_Ranked.png (Ranked comparison plot)
#   - results/Trait_Performance_Heatmap_Ranked.png (Trait vs Method performance heatmap)
#   - results/Method_Performance_Summary.csv (Comprehensive statistical summary table)
#
# Required R Packages: ggplot2, reshape2, PMCMRplus, dplyr, tidyr, multcompView, viridis
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


################################################################################
# STATISTICAL ANALYSIS: COMPARING TRAITS WITHIN EACH METHOD (ROBUST VERSION)
# Alternative approach with manual Nemenyi calculation
################################################################################

# Function to install and load packages
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    cat("Installing package:", package, "\n")
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
    cat("Package", package, "installed and loaded successfully!\n\n")
  } else {
    cat("Package", package, "already installed and loaded.\n")
  }
}

# List of required packages
required_packages <- c(
  "ggplot2",
  "reshape2",
  "dplyr",
  "tidyr",
  "ggpubr",
  "gridExtra",
  "knitr"
)

cat("================================================================================\n")
cat("CHECKING AND INSTALLING REQUIRED PACKAGES\n")
cat("================================================================================\n\n")

for (pkg in required_packages) {
  install_and_load(pkg)
}

cat("\n")
cat("================================================================================\n")
cat("   STATISTICAL ANALYSIS: COMPARING TRAITS WITHIN EACH METHOD\n")
cat("================================================================================\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC"), "\n")
cat("User:", Sys.info()["user"], "\n")
cat("================================================================================\n\n")

################################################################################
# FUNCTION TO CALCULATE KENDALL'S W
################################################################################

calculate_kendall_w <- function(data_matrix) {
  n <- nrow(data_matrix)  
  k <- ncol(data_matrix)  
  
  ranked_data <- t(apply(data_matrix, 1, rank))
  R_j <- colSums(ranked_data)
  R_bar <- mean(R_j)
  S <- sum((R_j - R_bar)^2)
  W <- (12 * S) / (n^2 * (k^3 - k))
  
  return(W)
}

################################################################################
# MANUAL NEMENYI POST-HOC TEST
################################################################################

nemenyi_test_manual <- function(data_matrix) {
  # data_matrix: rows = subjects, columns = treatments
  n <- nrow(data_matrix)  # number of subjects
  k <- ncol(data_matrix)  # number of treatments
  
  # Rank each row
  ranked_data <- t(apply(data_matrix, 1, rank))
  
  # Calculate mean ranks for each treatment
  mean_ranks <- colMeans(ranked_data)
  
  # Calculate pairwise differences
  trait_names <- colnames(data_matrix)
  p_matrix <- matrix(NA, nrow = k, ncol = k)
  rownames(p_matrix) <- trait_names
  colnames(p_matrix) <- trait_names
  
  # Critical value calculation for Nemenyi test
  # Standard error for rank differences
  SE <- sqrt((k * (k + 1)) / (6 * n))
  
  # Calculate p-values for all pairwise comparisons
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      # Absolute difference in mean ranks
      rank_diff <- abs(mean_ranks[i] - mean_ranks[j])
      
      # Calculate z-score
      z_score <- rank_diff / SE
      
      # Two-tailed p-value
      p_value <- 2 * (1 - pnorm(abs(z_score)))
      
      p_matrix[i, j] <- p_value
      p_matrix[j, i] <- p_value
    }
  }
  
  # Diagonal is 1 (same treatment)
  diag(p_matrix) <- 1
  
  return(list(
    p.value = p_matrix,
    mean.ranks = mean_ranks
  ))
}

################################################################################
# FUNCTION TO GENERATE COMPACT LETTER DISPLAY
################################################################################

generate_letter_groups <- function(p_matrix, trait_names, sig_level = 0.05) {
  n_traits <- length(trait_names)
  
  # Create adjacency matrix
  adj_matrix <- matrix(FALSE, n_traits, n_traits)
  rownames(adj_matrix) <- trait_names
  colnames(adj_matrix) <- trait_names
  
  for (i in 1:n_traits) {
    adj_matrix[i, i] <- TRUE
    for (j in 1:n_traits) {
      if (i != j && !is.na(p_matrix[i, j]) && p_matrix[i, j] >= sig_level) {
        adj_matrix[i, j] <- TRUE
        adj_matrix[j, i] <- TRUE
      }
    }
  }
  
  # Assign letters
  letters_assigned <- rep(NA, n_traits)
  names(letters_assigned) <- trait_names
  current_letter <- 1
  
  for (i in 1:n_traits) {
    if (is.na(letters_assigned[i])) {
      letters_assigned[i] <- current_letter
      
      for (j in (i+1):n_traits) {
        if (j <= n_traits && is.na(letters_assigned[j]) && adj_matrix[i, j]) {
          compatible <- TRUE
          for (k in 1:n_traits) {
            if (!is.na(letters_assigned[k]) && letters_assigned[k] == current_letter) {
              if (!adj_matrix[j, k]) {
                compatible <- FALSE
                break
              }
            }
          }
          if (compatible) {
            letters_assigned[j] <- current_letter
          }
        }
      }
      current_letter <- current_letter + 1
    }
  }
  
  letter_strings <- letters[letters_assigned]
  names(letter_strings) <- trait_names
  
  return(letter_strings)
}

################################################################################
# DATA STRUCTURE
################################################################################

cat("Inspecting data structure...\n")
cat("Column names in all_data:\n")
print(colnames(all_data))
cat("\n")

cat("First few rows of data:\n")
print(head(all_data))
cat("\n")

methods_vec <- unique(as.character(all_data$Method))
traits_vec <- unique(as.character(all_data$Trait))

cat("Methods analyzed:", paste(methods_vec, collapse = ", "), "\n")
cat("Traits compared:", paste(traits_vec, collapse = ", "), "\n\n")

id_col <- "CV_Iteration"
cat("Using '", id_col, "' as the replicate/fold identifier\n\n")

data_summary <- all_data %>%
  group_by(Method, Trait) %>%
  summarise(n_obs = n(), .groups = 'drop')

cat("Observations per Method-Trait combination:\n")
print(data_summary)
cat("\n")

################################################################################
# 1.FRIEDMAN TEST FOR EACH METHOD
################################################################################

cat("\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 1: FRIEDMAN TEST FOR EACH METHOD\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Null Hypothesis: No difference in predictive ability among traits within method\n")
cat("Alternative: At least one trait differs significantly within the method\n\n")

friedman_by_method <- list()
effect_sizes_by_method <- list()
posthoc_by_method <- list()
trait_rankings_by_method <- list()
letter_groups_by_method <- list()

friedman_summary_methods <- data.frame(
  Method = character(),
  Chi_Square = numeric(),
  df = integer(),
  p_value = character(),
  p_numeric = numeric(),
  Kendall_W = numeric(),
  Effect_Size = character(),
  Significant = character(),
  stringsAsFactors = FALSE
)

for (method in methods_vec) {
  cat("\n")
  cat("───────────────────────────────────────────────────────────────────────────\n")
  cat(paste0("METHOD: ", method, "\n"))
  cat("───────────────────────────────────────────────────────────────────────────\n")
  
  method_data <- all_data %>% 
    filter(Method == method) %>%
    select(all_of(c(id_col, "Trait", "Accuracy")))
  
  cat("  Data extracted:", nrow(method_data), "observations\n")
  
  method_wide <- method_data %>%
    pivot_wider(names_from = Trait, values_from = Accuracy, id_cols = all_of(id_col))
  
  cat("  Reshaped to wide format:", nrow(method_wide), "rows x", ncol(method_wide)-1, "trait columns\n")
  
  method_matrix <- method_wide %>%
    select(-all_of(id_col)) %>%
    as.matrix()
  
  if (any(is.na(method_matrix))) {
    cat("  WARNING: Missing values detected.Removing incomplete cases.\n")
    method_matrix <- method_matrix[complete.cases(method_matrix), ]
  }
  
  cat("  Final matrix for analysis:", nrow(method_matrix), "observations x", ncol(method_matrix), "traits\n")
  cat("  Trait columns:", paste(colnames(method_matrix), collapse = ", "), "\n")
  
  if (nrow(method_matrix) < 3) {
    cat("  ERROR: Insufficient data for Friedman test (need at least 3 observations)\n")
    next
  }
  
  # Friedman test
  friedman_result <- friedman.test(method_matrix)
  friedman_by_method[[method]] <- friedman_result
  
  # Kendall's W
  kendall_w <- calculate_kendall_w(method_matrix)
  effect_sizes_by_method[[method]] <- kendall_w
  
  effect_interpretation <- ifelse(kendall_w < 0.1, "Small",
                                  ifelse(kendall_w < 0.3, "Medium",
                                         ifelse(kendall_w < 0.5, "Large", "Very Large")))
  
  p_formatted <- ifelse(friedman_result$p.value < 0.001, "< 0.001",
                        ifelse(friedman_result$p.value < 0.01, sprintf("%.4f", friedman_result$p.value),
                               sprintf("%.3f", friedman_result$p.value)))
  
  sig_indicator <- ifelse(friedman_result$p.value < 0.001, "***",
                          ifelse(friedman_result$p.value < 0.01, "**",
                                 ifelse(friedman_result$p.value < 0.05, "*", "ns")))
  
  friedman_summary_methods <- rbind(friedman_summary_methods, data.frame(
    Method = method,
    Chi_Square = round(friedman_result$statistic, 3),
    df = friedman_result$parameter,
    p_value = p_formatted,
    p_numeric = friedman_result$p.value,
    Kendall_W = round(kendall_w, 4),
    Effect_Size = effect_interpretation,
    Significant = sig_indicator,
    stringsAsFactors = FALSE
  ))
  
  # Mean ranks
  trait_ranks <- colMeans(apply(method_matrix, 1, rank))
  trait_rankings_by_method[[method]] <- trait_ranks
  
  cat("  Friedman χ² =", round(friedman_result$statistic, 3), "\n")
  cat("  df =", friedman_result$parameter, "\n")
  cat("  p-value =", p_formatted, sig_indicator, "\n")
  cat("  Kendall's W =", round(kendall_w, 4), "(", effect_interpretation, "effect )\n")
  cat("  Mean Ranks (1 = best):\n")
  for (trait in names(sort(trait_ranks))) {
    cat("    ", trait, ":", round(trait_ranks[trait], 2), "\n")
  }
  
  # Nemenyi post-hoc test
  if (friedman_result$p.value < 0.05) {
    cat("\n  → Significant differences found.Performing Nemenyi post-hoc test...\n")
    
    # Use manual Nemenyi test
    posthoc_result <- nemenyi_test_manual(method_matrix)
    posthoc_by_method[[method]] <- posthoc_result
    
    p_matrix <- posthoc_result$p.value
    trait_names <- colnames(method_matrix)
    
    sig_level <- 0.05
    letter_vec <- generate_letter_groups(p_matrix, trait_names, sig_level)
    
    letter_groups_by_method[[method]] <- letter_vec
    
    cat("  Pairwise p-values (Nemenyi test):\n")
    print(round(p_matrix, 4))
    cat("\n  Significance groups (α = 0.05):\n")
    for (trait in names(sort(trait_ranks))) {
      cat("    ", trait, ":", letter_groups_by_method[[method]][trait], 
          " [Rank:", round(trait_ranks[trait], 2), "]\n")
    }
  } else {
    cat("\n  → No significant differences among traits (p ≥ 0.05)\n")
    letter_groups_by_method[[method]] <- setNames(rep("a", length(trait_ranks)), names(trait_ranks))
  }
}

################################################################################
# 2.SUMMARY TABLE
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 2: SUMMARY OF FRIEDMAN TESTS (All Methods)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

if (nrow(friedman_summary_methods) > 0) {
  print(kable(friedman_summary_methods[, -5], 
              col.names = c("Method", "χ²", "df", "p-value", "Kendall's W", 
                            "Effect Size", "Sig."),
              align = "lcccccc",
              format = "simple"))
  
  cat("\n")
  cat("Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant\n")
  cat("Kendall's W: < 0.1 (Small), 0.1-0.3 (Medium), 0.3-0.5 (Large), > 0.5 (Very Large)\n")
}

################################################################################
# 3.MEAN RANKS TABLE
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 3: MEAN RANKS OF TRAITS WITHIN EACH METHOD\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Lower rank = Better performance within that method\n\n")

if (length(trait_rankings_by_method) > 0) {
  trait_rank_table <- data.frame(Method = names(trait_rankings_by_method), stringsAsFactors = FALSE)
  
  for (trait in traits_vec) {
    ranks <- sapply(names(trait_rankings_by_method), function(m) {
      if (trait %in% names(trait_rankings_by_method[[m]])) {
        return(trait_rankings_by_method[[m]][trait])
      } else {
        return(NA)
      }
    })
    trait_rank_table[[trait]] <- round(ranks, 2)
  }
  
  print(kable(trait_rank_table,
              col.names = c("Method", traits_vec),
              align = paste0("l", paste(rep("c", length(traits_vec)), collapse = "")),
              format = "simple"))
}

################################################################################
# 4.SIGNIFICANCE GROUPS
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 4: SIGNIFICANCE GROUPS (Compact Letter Display)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Traits with same letter within a method are NOT significantly different (α = 0.05)\n\n")

if (length(letter_groups_by_method) > 0) {
  trait_letter_table <- data.frame(Method = names(letter_groups_by_method), stringsAsFactors = FALSE)
  
  for (trait in traits_vec) {
    letters_col <- sapply(names(letter_groups_by_method), function(m) {
      if (trait %in% names(letter_groups_by_method[[m]])) {
        return(letter_groups_by_method[[m]][trait])
      } else {
        return("")
      }
    })
    trait_letter_table[[trait]] <- letters_col
  }
  
  print(kable(trait_letter_table,
              col.names = c("Method", traits_vec),
              align = paste0("l", paste(rep("c", length(traits_vec)), collapse = "")),
              format = "simple"))
}

################################################################################
# 5.COMPREHENSIVE TABLE
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 5: COMPREHENSIVE PERFORMANCE SUMMARY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

if (length(trait_rankings_by_method) > 0 && length(letter_groups_by_method) > 0) {
  comprehensive_table <- data.frame(Method = names(trait_rankings_by_method), stringsAsFactors = FALSE)
  
  for (trait in traits_vec) {
    trait_summary <- all_data %>%
      filter(Trait == trait) %>%
      group_by(Method) %>%
      summarise(
        Mean = mean(Accuracy),
        SE = sd(Accuracy) / sqrt(n()),
        .groups = 'drop'
      )
    
    formatted_col <- sapply(names(trait_rankings_by_method), function(m) {
      mean_val <- trait_summary$Mean[trait_summary$Method == m]
      se_val <- trait_summary$SE[trait_summary$Method == m]
      
      if (trait %in% names(trait_rankings_by_method[[m]])) {
        rank_val <- trait_rankings_by_method[[m]][trait]
        letter_val <- letter_groups_by_method[[m]][trait]
        
        return(paste0(sprintf("%.3f", mean_val), " ± ", sprintf("%.3f", se_val),
                      " [", round(rank_val, 1), "] ", letter_val))
      } else {
        return(paste0(sprintf("%.3f", mean_val), " ± ", sprintf("%.3f", se_val)))
      }
    })
    
    comprehensive_table[[trait]] <- formatted_col
  }
  
  print(kable(comprehensive_table,
              col.names = c("Method", traits_vec),
              align = paste0("l", paste(rep("c", length(traits_vec)), collapse = "")),
              format = "simple"))
  
  cat("\n")
  cat("Format: Mean ± SE [Rank] Letter\n")
  cat("  - Mean: Average predictive ability\n")
  cat("  - SE: Standard error\n")
  cat("  - [Rank]: Mean rank (1 = best)\n")
  cat("  - Letter: Significance group\n")
}

################################################################################
# 6.DETAILED PAIRWISE COMPARISONS
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 6: DETAILED PAIRWISE COMPARISONS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

for (method in names(friedman_by_method)) {
  if (friedman_by_method[[method]]$p.value < 0.05) {
    cat("\n───────────────────────────────────────────────────────────────────────────\n")
    cat(paste0("METHOD: ", method, "\n"))
    cat("───────────────────────────────────────────────────────────────────────────\n\n")
    
    p_matrix <- posthoc_by_method[[method]]$p.value
    
    p_display <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
    rownames(p_display) <- rownames(p_matrix)
    colnames(p_display) <- colnames(p_matrix)
    
    for (i in 1:nrow(p_matrix)) {
      for (j in 1:ncol(p_matrix)) {
        if (is.na(p_matrix[i, j])) {
          p_display[i, j] <- "-"
        } else if (p_matrix[i, j] < 0.001) {
          p_display[i, j] <- "< 0.001 ***"
        } else if (p_matrix[i, j] < 0.01) {
          p_display[i, j] <- sprintf("%.4f **", p_matrix[i, j])
        } else if (p_matrix[i, j] < 0.05) {
          p_display[i, j] <- sprintf("%.4f *", p_matrix[i, j])
        } else {
          p_display[i, j] <- sprintf("%.4f ns", p_matrix[i, j])
        }
      }
    }
    
    print(kable(p_display, align = "c", format = "simple"))
    cat("\n")
  }
}

################################################################################
# 7.EXPORT RESULTS
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 7: EXPORTING RESULTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

if (nrow(friedman_summary_methods) > 0) {
  write.csv(friedman_summary_methods, 
            paste0(dir, "/Friedman_Test_Traits_Within_Methods.csv"), 
            row.names = FALSE)
  cat("✓ Saved: Friedman_Test_Traits_Within_Methods.csv\n")
  
  if (exists("trait_rank_table")) {
    write.csv(trait_rank_table, 
              paste0(dir, "/Trait_Ranks_Within_Methods.csv"), 
              row.names = FALSE)
    cat("✓ Saved: Trait_Ranks_Within_Methods.csv\n")
  }
  
  if (exists("trait_letter_table")) {
    write.csv(trait_letter_table, 
              paste0(dir, "/Significance_Letters_Traits_Within_Methods.csv"), 
              row.names = FALSE)
    cat("✓ Saved: Significance_Letters_Traits_Within_Methods.csv\n")
  }
  
  if (exists("comprehensive_table")) {
    write.csv(comprehensive_table, 
              paste0(dir, "/Comprehensive_Traits_Within_Methods.csv"), 
              row.names = FALSE)
    cat("✓ Saved: Comprehensive_Traits_Within_Methods.csv\n")
  }
  
  for (method in names(posthoc_by_method)) {
    p_matrix <- posthoc_by_method[[method]]$p.value
    write.csv(p_matrix, 
              paste0(dir, "/Nemenyi_Pairwise_", gsub(" ", "_", method), "_Traits.csv"))
    cat(paste0("✓ Saved: Nemenyi_Pairwise_", gsub(" ", "_", method), "_Traits.csv\n"))
  }
}

################################################################################
# 8.INTERPRETATION
################################################################################

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("STEP 8: KEY FINDINGS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

if (nrow(friedman_summary_methods) > 0) {
  n_sig <- sum(friedman_summary_methods$p_numeric < 0.05)
  cat("SUMMARY:\n")
  cat("  • Methods with significant trait differences:", n_sig, "out of", nrow(friedman_summary_methods), "\n")
  cat("  • Methods without significant differences:", nrow(friedman_summary_methods) - n_sig, "\n\n")
  
  if (n_sig > 0) {
    cat("METHODS WITH SIGNIFICANT DIFFERENCES (p < 0.05):\n")
    sig_methods <- friedman_summary_methods$Method[friedman_summary_methods$p_numeric < 0.05]
    for (m in sig_methods) {
      cat("  •", m, "(p", 
          ifelse(friedman_summary_methods$p_numeric[friedman_summary_methods$Method == m] < 0.001,
                 "< 0.001", 
                 paste("=", sprintf("%.4f", friedman_summary_methods$p_numeric[friedman_summary_methods$Method == m]))),
          ", W =", sprintf("%.3f", friedman_summary_methods$Kendall_W[friedman_summary_methods$Method == m]),
          ")\n")
    }
  }
}

cat("\n")
cat("================================================================================\n")
cat("                    ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
################################################################################
# VISUALIZATION: TRAIT PERFORMANCE WITHIN EACH METHOD
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

cat("\n")
cat("================================================================================\n")
cat("CREATING VISUALIZATION OF TRAIT PERFORMANCE WITHIN METHODS\n")
cat("================================================================================\n\n")

################################################################################
# PREPARE DATA FOR PLOTTING
################################################################################

# Extract mean, SE, ranks, and letters from comprehensive_table
plot_data <- data.frame(
  Method = character(),
  Trait = character(),
  Mean = numeric(),
  SE = numeric(),
  Rank = numeric(),
  Letter = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(comprehensive_table)) {
  method <- comprehensive_table$Method[i]
  
  for (trait in traits_vec) {
    # Parse the formatted string: "Mean ± SE [Rank] Letter"
    value_str <- comprehensive_table[[trait]][i]
    
    # Extract mean
    mean_val <- as.numeric(gsub(" ±.*", "", value_str))
    
    # Extract SE
    se_val <- as.numeric(gsub(".*± ([0-9.]+).*", "\\1", value_str))
    
    # Extract rank (if present)
    if (grepl("\\[", value_str)) {
      rank_val <- as.numeric(gsub(".*\\[([0-9.]+)\\].*", "\\1", value_str))
      letter_val <- gsub(".*\\] ", "", value_str)
    } else {
      rank_val <- NA
      letter_val <- ""
    }
    
    plot_data <- rbind(plot_data, data.frame(
      Method = method,
      Trait = trait,
      Mean = mean_val,
      SE = se_val,
      Rank = rank_val,
      Letter = letter_val,
      stringsAsFactors = FALSE
    ))
  }
}

# Reorder methods by overall performance (average across traits)
method_order <- plot_data %>%
  group_by(Method) %>%
  summarise(avg_mean = mean(Mean), .groups = 'drop') %>%
  arrange(desc(avg_mean)) %>%
  pull(Method)

plot_data$Method <- factor(plot_data$Method, levels = method_order)

# Rename traits for display
plot_data$Trait_Display <- factor(plot_data$Trait,
                                  levels = c("Is_Dead", "DPC_date", "DPC_time"),
                                  labels = c("Binary Survival", "DPC_Date", "DPC_Time"))

# Get Kendall's W values for legend
kendall_w_values <- friedman_summary_methods %>%
  select(Method, Kendall_W) %>%
  arrange(match(Method, method_order))

cat("Data prepared for plotting:\n")
print(head(plot_data))
cat("\n")

################################################################################
# CREATE MAIN PLOT WITH STATISTICAL ANNOTATIONS
################################################################################

create_trait_comparison_plot <- function() {
  
  # Color palette
  trait_colors <- c(
    "Binary Survival" = "#E31A1C",  # Red
    "DPC_Date" = "#1F78B4",          # Blue
    "DPC_Time" = "#33A02C"           # Green
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Method, y = Mean, fill = Trait_Display)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85), 
             color = "black", alpha = 0.85, width = 0.75, linewidth = 0.8) +
    
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(width = 0.85), 
                  width = 0.2, linewidth = 0.6) +
    
    # Add rank numbers in RED inside bars
    geom_text(aes(label = ifelse(!is.na(Rank), sprintf("%.1f", Rank), ""), 
                  y = Mean/2),
              position = position_dodge(width = 0.85), 
              color = "red", size = 5, fontface = "bold") +
    
    # Add significance letters in RED above error bars
    geom_text(aes(label = Letter, y = Mean + SE + 0.015),
              position = position_dodge(width = 0.85), 
              vjust = 0, size = 5.5, fontface = "bold", 
              color = "red") +
    
    scale_fill_manual(
      values = trait_colors,
      name = "Trait",
      labels = c(
        paste0("Binary Survival (W = ", sprintf("%.3f", mean(kendall_w_values$Kendall_W)), ")"),
        paste0("DPC_Date (W = ", sprintf("%.3f", mean(kendall_w_values$Kendall_W)), ")"),
        paste0("DPC_Time (W = ", sprintf("%.3f", mean(kendall_w_values$Kendall_W)), ")")
      )
    ) +
    
    labs(
      title = "Mean Predictive Ability by Method with Statistical Significance",
      subtitle = "Red letters: significance groups within each method (p < 0.05) | Red numbers: mean ranks (1 = best)\nFriedman test followed by Nemenyi post-hoc test",
      x = "Prediction Method",
      y = "Mean Predictive Ability ± SE",
      caption = "Traits with the same letter within a method are not significantly different"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15), 
                                   color = "darkred", face = "bold"),
      plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", margin = margin(t = 10)),
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.key.size = unit(1.2, "cm"),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray92", linewidth = 0.3),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    
    scale_y_continuous(breaks = seq(0, 0.6, 0.05), 
                       limits = c(0, max(plot_data$Mean + plot_data$SE) + 0.05),
                       labels = scales::number_format(accuracy = 0.01))
  
  return(p)
}

################################################################################
# CREATE ALTERNATIVE HEATMAP VISUALIZATION
################################################################################

create_heatmap_plot <- function() {
  
  # Prepare data for heatmap
  heatmap_data <- plot_data %>%
    select(Method, Trait, Mean, Letter, Rank)
  
  p_heat <- ggplot(heatmap_data, aes(x = Trait, y = Method, fill = Mean)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%.3f\n[%.1f] %s", Mean, Rank, Letter)),
              size = 4, fontface = "bold", color = "black") +
    
    scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850",
                         midpoint = median(heatmap_data$Mean),
                         name = "Mean\nAccuracy") +
    
    labs(
      title = "Heatmap: Trait Performance Within Each Method",
      subtitle = "Values show: Mean [Rank] Letter | Letter = significance group",
      x = "Trait",
      y = "Prediction Method"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 15)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()
    )
  
  return(p_heat)
}

################################################################################
# CREATE FACETED PLOT BY TRAIT
################################################################################

create_faceted_plot <- function() {
  
  p_facet <- ggplot(plot_data, aes(x = Method, y = Mean, fill = Trait_Display)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.85, width = 0.7) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  width = 0.3, linewidth = 0.6) +
    
    # Rank in red
    geom_text(aes(label = sprintf("%.1f", Rank), y = Mean/2),
              color = "red", size = 4, fontface = "bold") +
    
    # Letter in red
    geom_text(aes(label = Letter, y = Mean + SE + 0.01),
              color = "red", size = 5, fontface = "bold", vjust = 0) +
    
    facet_wrap(~Trait_Display, ncol = 1, scales = "free_y") +
    
    scale_fill_manual(values = c(
      "Binary Survival" = "#E31A1C",
      "DPC_Date" = "#1F78B4",
      "DPC_Time" = "#33A02C"
    )) +
    
    labs(
      title = "Trait Performance Across Methods (Faceted View)",
      subtitle = "Red numbers: ranks | Red letters: significance groups",
      x = "Prediction Method",
      y = "Mean Predictive Ability ± SE"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkred", face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  return(p_facet)
}

################################################################################
# GENERATE AND SAVE PLOTS
################################################################################

cat("Creating visualizations...\n\n")

# Main grouped bar chart
main_plot <- create_trait_comparison_plot()
print(main_plot)

ggsave(paste0(dir, "/Trait_Comparison_Within_Methods_Grouped_Bars.png"), 
       main_plot, width = 16, height = 10, dpi = 300)
cat("✓ Saved: Trait_Comparison_Within_Methods_Grouped_Bars.png\n")

# Heatmap
heatmap_plot <- create_heatmap_plot()
print(heatmap_plot)

ggsave(paste0(dir, "/Trait_Comparison_Within_Methods_Heatmap.png"), 
       heatmap_plot, width = 10, height = 12, dpi = 300)
cat("✓ Saved: Trait_Comparison_Within_Methods_Heatmap.png\n")

# Faceted plot
faceted_plot <- create_faceted_plot()
print(faceted_plot)

ggsave(paste0(dir, "/Trait_Comparison_Within_Methods_Faceted.png"), 
       faceted_plot, width = 14, height = 12, dpi = 300)
cat("✓ Saved: Trait_Comparison_Within_Methods_Faceted.png\n")

################################################################################
# PUBLICATION-READY VERSION (HIGH RESOLUTION)
################################################################################

ggsave(paste0(dir, "/Publication_Trait_Comparison_Within_Methods.png"), 
       main_plot, width = 18, height = 11, dpi = 600)
cat("✓ Saved: Publication_Trait_Comparison_Within_Methods.png (high resolution)\n")

ggsave(paste0(dir, "/Publication_Trait_Comparison_Within_Methods.pdf"), 
       main_plot, width = 18, height = 11, device = "pdf")
cat("✓ Saved: Publication_Trait_Comparison_Within_Methods.pdf (vector format)\n")

cat("\n")
cat("================================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Three types of plots created:\n")
cat("1.Grouped Bar Chart - Shows all traits side-by-side for each method\n")
cat("2.Heatmap - Shows performance matrix with color coding\n")
cat("3.Faceted Plot - Separate panels for each trait\n\n")

cat("Key features:\n")
cat("  • Red numbers inside bars = Mean ranks (1 = best trait for that method)\n")
cat("  • Red letters above bars = Significance groups (same letter = not significantly different)\n")
cat("  • Error bars = Standard error\n")
cat("  • Methods ordered by overall performance (left = best)\n\n")
