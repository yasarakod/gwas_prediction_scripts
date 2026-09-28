# ============================================================================== #
# Script: 15_population_size_friedman_analysis.R
# Title: Training Population Size Performance & Friedman Testing
# Description:
#   Evaluates prediction metrics across training population sizes (50 to 474 individuals) and CV folds (3-fold vs 5-fold), generating Friedman test comparisons and boxplots.
#
# Inputs:
#   - results/GBLUP_Folds_*_Times_10_* (Cross-validation output files for different sample sizes)
#
# Outputs:
#   - results/Table_<pmethod><namef>.csv (Summary table)
#   - results/Friedman_Results_<pmethod><namef>.csv (Friedman test results)
#   - results/PopFolds_Rankings_<pmethod><namef>.csv (Population size rankings)
#   - results/Accuracy_Plot2_<pmethod><namef>.png (Accuracy multi-panel plot)
#   - results/Ability_Plot2_<pmethod><namef>.png (Predictive ability multi-panel plot)
#
# Required R Packages: ggplot2, reshape2, PMCMRplus, dplyr
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
library(ggplot2)
library(reshape2)
library(PMCMRplus)
library(dplyr)

# Function to perform Friedman test and generate ranking letters
perform_friedman_analysis <- function(data_matrix) {
  
  # Friedman test
  friedman_result <- friedman.test(data_matrix)
  cat("=== FRIEDMAN TEST RESULTS ===\n")
  print(friedman_result)
  
  # Calculate effect size (Kendall's W)
  n_blocks <- nrow(data_matrix)
  n_treatments <- ncol(data_matrix)
  chi_squared <- friedman_result$statistic
  kendall_w <- chi_squared / (n_blocks * (n_treatments - 1))
  
  cat("\nEffect size (Kendall's W):", round(kendall_w, 4), "\n")
  
  # Initialize ranking letters
  ranking_letters <- rep("a", ncol(data_matrix))
  names(ranking_letters) <- colnames(data_matrix)
  
  # If Friedman test is significant, perform post-hoc analysis
  if(friedman_result$p.value < 0.05) {
    cat("\nFriedman test is significant (p < 0.05). Performing post-hoc analysis...\n")
    
    # Nemenyi post-hoc test
    nemenyi_result <- tryCatch({
      PMCMRplus::frdAllPairsNemenyiTest(data_matrix)
    }, error = function(e) {
      cat("Nemenyi test failed:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(nemenyi_result)) {
      cat("Nemenyi test completed successfully.\n")
      
      # Calculate mean ranks (lower rank = better performance since higher values are better)
      ranks <- apply(data_matrix, 1, function(x) rank(-x))  # Negative for descending rank
      mean_ranks <- rowMeans(ranks)
      
      # Sort methods by mean rank
      sorted_methods <- names(sort(mean_ranks))
      
      # Generate ranking letters based on significant differences
      p_matrix <- nemenyi_result$p.value
      letters <- generate_ranking_letters_fixed(p_matrix, mean_ranks)
      ranking_letters <- letters[colnames(data_matrix)]
      
    } else {
      cat("Post-hoc analysis failed.Using simple ranking based on means.\n")
      # Fallback: simple ranking based on means
      mean_values <- colMeans(data_matrix)
      sorted_methods <- names(sort(mean_values, decreasing = TRUE))
      
      # Assign letters based on rank (best gets 'a')
      for(i in 1:length(sorted_methods)) {
        ranking_letters[sorted_methods[i]] <- letters[i]
      }
    }
    
  } else {
    cat("\nFriedman test is not significant (p ≥ 0.05).\n")
    cat("All methods have similar performance - assigning same letter.\n")
    # All methods get the same letter if no significant difference
    ranking_letters <- rep("a", length(ranking_letters))
    names(ranking_letters) <- colnames(data_matrix)
  }
  
  return(list(
    friedman_result = friedman_result,
    kendall_w = kendall_w,
    ranking_letters = ranking_letters
  ))
}

# Fixed function to generate ranking letters from p-value matrix
generate_ranking_letters_fixed <- function(p_matrix, mean_ranks, alpha = 0.05) {
  methods <- names(mean_ranks)
  n_methods <- length(methods)
  
  # Check if p_matrix is valid
  if(is.null(p_matrix) || !is.matrix(p_matrix) || n_methods <= 1) {
    cat("Invalid p-value matrix or insufficient methods.Using mean-based ranking.\n")
    return(assign_letters_by_mean(mean_ranks))
  }
  
  # Create adjacency matrix for significant differences
  sig_diff <- matrix(FALSE, n_methods, n_methods)
  rownames(sig_diff) <- colnames(sig_diff) <- methods
  
  # Fill significance matrix
  for(i in 1:(n_methods-1)) {
    for(j in (i+1):n_methods) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      p_val <- tryCatch({
        p_matrix[method1, method2]
      }, error = function(e) {
        tryCatch({
          p_matrix[method2, method1]
        }, error = function(e2) {
          return(NA)
        })
      })
      
      if(!is.na(p_val) && p_val < alpha) {
        sig_diff[method1, method2] <- TRUE
        sig_diff[method2, method1] <- TRUE
      }
    }
  }
  
  # Sort methods by mean rank for letter assignment
  sorted_methods <- names(sort(mean_ranks))
  
  # Initialize letters
  letters_assigned <- rep(NA, n_methods)
  names(letters_assigned) <- methods
  current_letter <- 1
  
  for(method in sorted_methods) {
    if(is.na(letters_assigned[method])) {
      # Find all methods not significantly different from current method
      group_methods <- c(method)
      
      for(other_method in sorted_methods) {
        if(other_method != method && is.na(letters_assigned[other_method])) {
          if(!sig_diff[method, other_method]) {
            group_methods <- c(group_methods, other_method)
          }
        }
      }
      
      # Assign same letter to all methods in this group
      if(current_letter <= 26) {
        letter <- letters[current_letter]
      } else {
        letter <- paste0(letters[(current_letter-1) %% 26 + 1], ceiling(current_letter/26))
      }
      
      for(m in group_methods) {
        letters_assigned[m] <- letter
      }
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# Fallback function to assign letters based on means only
assign_letters_by_mean <- function(mean_ranks) {
  methods <- names(mean_ranks)
  sorted_methods <- names(sort(mean_ranks))
  
  letters_assigned <- rep(NA, length(methods))
  names(letters_assigned) <- methods
  
  # Simple ranking: each method gets a different letter based on rank
  for(i in 1:length(sorted_methods)) {
    if(i <= 26) {
      letter <- letters[i]
    } else {
      letter <- paste0(letters[(i-1) %% 26 + 1], ceiling(i/26))
    }
    letters_assigned[sorted_methods[i]] <- letter
  }
  
  return(letters_assigned)
}

# Parse Arguments
dir <- "data"
graph_type <- "D"  # New graph type for folds comparison
trait_name <- c("Is_Dead")
pmethods <- c("GBLUP")
iter <- 10
pop_sizes <- c("50","100", "150", "200", "250", "300", "350", "400", "474")
marker_sizes <- c("top_1000")
Fixed_effect <- c("NULL")
h2 <- 1

# Define fold comparisons
folds_comparison <- c(3, 5)

data_path = paste0(dir, "/")

if (graph_type=="D"){
  
  sname <- trait_name
  pmethod <- pmethods
  
  # Create comparison matrix for all combinations of population sizes and folds
  total_configs <- length(pop_sizes) * length(folds_comparison)
  compareM <- matrix(0, iter, total_configs)
  
  # Create configuration names
  config_names <- c()
  config_counter <- 1
  
  # Process each combination of population size and fold
  for (pop_idx in seq_along(pop_sizes)) {
    for (fold_idx in seq_along(folds_comparison)) {
      pop_size <- pop_sizes[pop_idx]
      fold_val <- folds_comparison[fold_idx]
      
      filename <- paste0(data_path, pmethod, "_Folds_", fold_val, "_Times_", iter, "_", sname, "_", Fixed_effect, "_", pop_size, "_", marker_sizes, ".txt")
      
      cat("Looking for file:", filename, "\n")
      
      # Check if file exists
      if (file.exists(filename)) {
        cat("File found, reading data...\n")
        datap <- read.table(filename)
        for (n in 1:iter) {
          compareM[n, config_counter] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
        }
      } else {
        cat("Warning: File not found:", filename, "\n")
        compareM[, config_counter] <- NA  # Fill with NA if file doesn't exist
      }
      
      # Create configuration name
      config_names <- c(config_names, paste0("Pop", pop_size, "_", fold_val, "Fold"))
      config_counter <- config_counter + 1
    }
  }
  
  # Set column names and row names
  colnames(compareM) <- config_names
  rownames(compareM) <- paste0("CV", 1:iter)
  
  # Remove columns with all NA values
  compareM <- compareM[, !apply(is.na(compareM), 2, all)]
  
  # Check if we have valid data
  if (ncol(compareM) == 0) {
    stop("No valid data files found.Please check file paths and names.")
  }
  
  namef <- paste0("_Times_", iter, "_", sname, "_", Fixed_effect, "_", marker_sizes, "_pop_folds_comparison")
  
  # Save the comparison matrix
  write.table(compareM, paste0(data_path,"Table_",pmethod,namef,".csv"), sep = ",", quote = FALSE)
  
  # Perform statistical analysis if we have more than one column
  if (ncol(compareM) > 1) {
    # Perform Friedman test analysis
    friedman_analysis <- perform_friedman_analysis(compareM)
    ranking_letters <- friedman_analysis$ranking_letters
    
    # Save Friedman test results
    friedman_summary <- data.frame(
      Test_Statistic = friedman_analysis$friedman_result$statistic,
      P_Value = friedman_analysis$friedman_result$p.value,
      Kendall_W = friedman_analysis$kendall_w,
      Significant = friedman_analysis$friedman_result$p.value < 0.05
    )
    write.csv(friedman_summary, paste0(data_path, "Friedman_Results_", pmethod, namef, ".csv"), row.names = FALSE)
    
    # Create ranking summary table
    ranking_summary <- data.frame(
      Configuration = names(ranking_letters),
      Mean_Prediction_Ability = round(colMeans(compareM, na.rm = TRUE), 4),
      SD = round(apply(compareM, 2, sd, na.rm = TRUE), 4),
      Ranking_Letter = ranking_letters[names(ranking_letters)]
    ) %>%
      arrange(desc(Mean_Prediction_Ability))
    
    write.csv(ranking_summary, paste0(data_path, "PopFolds_Rankings_", pmethod, namef, ".csv"), row.names = FALSE)
  }
  
  # Melt the data for plotting
  plot_data <- reshape2::melt(compareM)
  
  # Extract population size and fold information
  plot_data$pop_size <- as.numeric(gsub("Pop([0-9]+)_.*", "\\1", plot_data$Var2))
  plot_data$fold_number <- as.numeric(gsub(".*_([0-9]+)Fold", "\\1", plot_data$Var2))
  plot_data$group <- paste0(plot_data$fold_number, "-Fold CV")
  
  # Create proper ordering for x-axis (by population size)
  plot_data$x_label <- factor(plot_data$pop_size, levels = sort(unique(plot_data$pop_size)))
  
  # Calculate y-positions for ranking letters
  max_values <- plot_data %>%
    group_by(x_label, group) %>%
    summarise(max_val = max(value, na.rm = TRUE), .groups = 'drop')
  
  # Add ranking letters to max_values
  if (exists("ranking_letters")) {
    max_values$ranking_letter <- ""
    for (i in 1:nrow(max_values)) {
      pop_val <- as.character(max_values$x_label[i])
      fold_val <- gsub("-Fold CV", "", max_values$group[i])
      config_name <- paste0("Pop", pop_val, "_", fold_val, "Fold")
      
      if (config_name %in% names(ranking_letters)) {
        max_values$ranking_letter[i] <- ranking_letters[config_name]
      }
    }
  } else {
    max_values$ranking_letter <- ""
  }
  
  # Create the first plot - Prediction Accuracy (with h2 correction)
  # Create the first plot - Prediction Accuracy (with h2 correction) - WITH CONNECTED MEANS
  p1 <- ggplot(plot_data, aes(x=x_label, y=value/sqrt(h2), fill=group)) + 
    geom_boxplot(position = position_dodge(width = 0.8), color="black", alpha = 0.8) +
    scale_fill_manual(
      values = c("3-Fold CV" = "#E31A1C", "5-Fold CV" = "#1F78B4"),
      name = "CV Strategy"
    ) +
    labs(title=paste0("Prediction Accuracy: Population Size vs CV Folds - ", sname), 
         x="Population Size", 
         y="Prediction Accuracy") +
    # Add mean points
    stat_summary(aes(group = group), fun=mean, geom="point", shape=21, size=2, 
                 color="black", fill="yellow", stroke=1,
                 position = position_dodge(width = 0.8)) +
    # Add dashed lines connecting mean points for each CV strategy
    stat_summary(aes(group = group, color = group), fun=mean, geom="line", 
                 linetype="dashed", size=1.2, 
                 position = position_dodge(width = 0.8)) +
    # Set line colors to match the fill colors
    scale_color_manual(
      values = c("3-Fold CV" = "#E31A1C", "5-Fold CV" = "#1F78B4"),
      guide = "none"  # Hide the color legend since we already have fill legend
    ) +
    theme(
      legend.title = element_text(colour="blue", size=16, face="bold", margin=margin(b=15)),
      legend.text = element_text(size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=18, face="bold", margin=margin(b=30)),
      axis.title.x = element_text(size=18, face="bold", margin=margin(t=20)),
      axis.title.y = element_text(size=18, face="bold", margin=margin(r=20)),
      axis.text = element_text(size=12, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1, size=12, face="bold", margin=margin(t=10)),
      axis.text.y = element_text(size=12, face="bold", margin=margin(r=10)),
      plot.background = element_rect(fill = "#FFFFFF"),
      legend.background = element_rect(linetype = "dashed"),
      panel.border = element_rect(fill = NA),
      plot.margin = margin(t=40, r=20, b=30, l=30)
    )
  
  ggsave(paste0(data_path,"Accuracy_Plot2_",pmethod,namef,".png"), p1, scale=1, dpi = 300, height = 8, width = 14)
  
  # Create the second plot - Predictive Ability (without h2 correction) - WITH CONNECTED MEANS
  p2 <- ggplot(plot_data, aes(x=x_label, y=value, fill=group)) + 
    geom_boxplot(position = position_dodge(width = 0.8), color="black", alpha = 0.8) +
    scale_fill_manual(
      values = c("3-Fold CV" = "#E31A1C", "5-Fold CV" = "#1F78B4"),
      name = "CV Strategy"
    ) +
    labs(title=paste0("GBLUP - Binary Survival"), 
         x="Population Size", 
         y="Predictive Ability") +
    # Add mean points
    stat_summary(aes(group = group), fun=mean, geom="point", shape=21, size=2, 
                 color="black", fill="yellow", stroke=1,
                 position = position_dodge(width = 0.8)) +
    # Add dashed lines connecting mean points for each CV strategy
    stat_summary(aes(group = group, color = group), fun=mean, geom="line", 
                 linetype="dashed", size=1.2, 
                 position = position_dodge(width = 0.8)) +
    # Set line colors to match the fill colors
    scale_color_manual(
      values = c("3-Fold CV" = "#E31A1C", "5-Fold CV" = "#1F78B4"),
      guide = "none"  # Hide the color legend since we already have fill legend
    ) +
    theme(
      legend.title = element_text(colour="blue", size=16, face="bold", margin=margin(b=15)),
      legend.text = element_text(size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=18, face="bold", margin=margin(b=30)),
      axis.title.x = element_text(size=18, face="bold", margin=margin(t=20)),
      axis.title.y = element_text(size=18, face="bold", margin=margin(r=20)),
      axis.text = element_text(size=12, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1, size=16, face="bold", margin=margin(t=10)),
      axis.text.y = element_text(size=16, face="bold", margin=margin(r=10)),
      plot.background = element_rect(fill = "#FFFFFF"),
      legend.background = element_rect(linetype = "dashed"),
      panel.border = element_rect(fill = NA),
      plot.margin = margin(t=40, r=20, b=30, l=30)
    )
  
  ggsave(paste0(data_path,"Ability_Plot2_",pmethod,namef,".png"), p2, scale=1, dpi = 300, height = 8, width = 10)
  
  
  # Print summary to console
  cat("\n=== POPULATION SIZE vs CV FOLDS COMPARISON ANALYSIS COMPLETE ===\n")
  cat("Files generated:\n")
  cat("- Accuracy_Plot_with_Rankings_", pmethod, namef, ".png\n")
  cat("- Ability_Plot_with_Rankings_", pmethod, namef, ".png\n")
  
  if (exists("ranking_letters") && length(ranking_letters) > 1) {
    cat("- Friedman_Results_", pmethod, namef, ".csv\n")
    cat("- PopFolds_Rankings_", pmethod, namef, ".csv\n")
    
    cat("\nConfiguration Ranking Letters:\n")
    for(i in 1:length(ranking_letters)) {
      cat(names(ranking_letters)[i], ": ", ranking_letters[i], "\n")
    }
    
    cat("\nNote: Configurations sharing the same letter are not significantly different.\n")
    cat("Letters are assigned alphabetically with 'a' being the best performing group.\n")
  } else {
    cat("Note: Insufficient configurations for statistical comparison.\n")
  }
  
  cat("\nExpected file pattern:\n")
  cat("GBLUP_Folds_[3|5]_Times_10_DPC_time_NULL_[pop_size]_top_1000.txt\n")
}
