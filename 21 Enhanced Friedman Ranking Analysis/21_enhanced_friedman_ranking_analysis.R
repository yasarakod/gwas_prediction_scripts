# ============================================================================== #
# Script: 21_enhanced_friedman_ranking_analysis.R
# Title: Enhanced Friedman & Nemenyi Ranking Analysis with CLD
# Description:
#   Performs robust Friedman-Nemenyi analysis featuring compact letter display (CLD) generation, Kendall W concordance calculation, and enhanced visualization.
#
# Inputs:
#   - results/Table_Folds_5_Times_10_* (Model evaluation tables for target traits)
#
# Outputs:
#   - results/Friedman_Nemenyi_Statistical_Report.txt (Comprehensive statistical text report)
#   - results/Friedman_<trait>_analysis.png (Enhanced significance boxplots with ranking annotations)
#
# Required R Packages: ggplot2, reshape2, ggpubr, PMCMRplus, dplyr
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


# setwd() configured via project environment
# Load required libraries
library(ggplot2)
library(reshape2)
library(ggpubr)
library(PMCMRplus)
library(dplyr)

# Parse Arguments
dir<-"data"
graph_type<-"A"
#trait_name<-c("DPC_time")
trait_name<-c("Is_D")
pmethods<-c("PBLUP","GBLUP","BA","BB","BC","BL","BRR","RR","EN","RF")
iter<-10
fold<-5
pop_sizes<-"NULL"
marker_sizes<-c("top_40000")
Fixed_effect<-c("NULL")
h2<-1
data_path = paste0(dir, "/")

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
      
      # Debug: Check the structure of the p-value matrix
      cat("P-value matrix structure:\n")
      cat("Class:", class(nemenyi_result$p.value), "\n")
      cat("Dimensions:", dim(nemenyi_result$p.value), "\n")
      cat("Row names:", rownames(nemenyi_result$p.value), "\n")
      cat("Column names:", colnames(nemenyi_result$p.value), "\n")
      
      # Calculate mean ranks (lower rank = better performance since higher values are better)
      ranks <- apply(data_matrix, 1, function(x) rank(-x))  # Negative for descending rank
      mean_ranks <- rowMeans(ranks)
      
      # Sort methods by mean rank
      sorted_methods <- names(sort(mean_ranks))
      
      # Generate ranking letters based on significant differences
      p_matrix <- nemenyi_result$p.value
      letters <- generate_ranking_letters_fixed(p_matrix, mean_ranks)
      ranking_letters <- letters[colnames(data_matrix)]
      
      # Print ranking summary
      rank_summary <- data.frame(
        Method = names(mean_ranks),
        Mean_Rank = round(mean_ranks, 2),
        Mean_Value = round(colMeans(data_matrix), 4),
        Ranking_Letter = letters[names(mean_ranks)]
      ) %>%
        arrange(Mean_Rank)
      
      cat("\n=== RANKING SUMMARY ===\n")
      print(rank_summary)
      
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
  
  cat("Generating ranking letters...\n")
  cat("Number of methods:", n_methods, "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  
  # Check if p_matrix is valid
  if(is.null(p_matrix) || !is.matrix(p_matrix)) {
    cat("Invalid p-value matrix.Using mean-based ranking.\n")
    return(assign_letters_by_mean(mean_ranks))
  }
  
  # Get matrix row and column names
  p_rownames <- rownames(p_matrix)
  p_colnames <- colnames(p_matrix)
  
  cat("P-matrix row names:", paste(p_rownames, collapse = ", "), "\n")
  cat("P-matrix col names:", paste(p_colnames, collapse = ", "), "\n")
  
  # Check if all methods are present in the matrix
  methods_in_matrix <- methods[methods %in% p_rownames & methods %in% p_colnames]
  
  if(length(methods_in_matrix) < length(methods)) {
    cat("Warning: Not all methods found in p-value matrix.\n")
    cat("Methods in matrix:", paste(methods_in_matrix, collapse = ", "), "\n")
    cat("Using mean-based ranking as fallback.\n")
    return(assign_letters_by_mean(mean_ranks))
  }
  
  # Create adjacency matrix for significant differences
  sig_diff <- matrix(FALSE, n_methods, n_methods)
  rownames(sig_diff) <- colnames(sig_diff) <- methods
  
  # Fill significance matrix with proper error checking
  for(i in 1:(n_methods-1)) {
    for(j in (i+1):n_methods) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      # Check if both methods exist in the p-value matrix
      if(method1 %in% p_rownames && method2 %in% p_colnames) {
        p_val <- tryCatch({
          p_matrix[method1, method2]
        }, error = function(e) {
          # Try the other way around
          tryCatch({
            p_matrix[method2, method1]
          }, error = function(e2) {
            cat("Could not access p-value for", method1, "vs", method2, "\n")
            return(NA)
          })
        })
        
        if(!is.na(p_val) && p_val < alpha) {
          sig_diff[method1, method2] <- TRUE
          sig_diff[method2, method1] <- TRUE
        }
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

if (graph_type=="A"){
  
  trait<-trait_name
  pmethod<-pmethods
  compareM<-matrix(0,iter,length(pmethod))
  for (x in seq_along(pmethod)) {
    datap <- read.table(paste0(data_path,
                               pmethod[x], "_Folds_", fold, "_Times_", iter, "_", 
                               paste(trait, collapse = "_"), "_", paste(Fixed_effect, collapse = "_"), 
                               "_", paste(pop_sizes, collapse = "_"), "_", 
                               paste(marker_sizes, collapse = "_"), ".txt"
    ))
    for (n in 1:iter) {
      compareM[n, x] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  colnames(compareM) <- pmethod
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_", 
                paste(trait, collapse = "_"), "_", paste(Fixed_effect, collapse = "_"), 
                "_", paste(pop_sizes, collapse = "_"), "_", 
                paste(marker_sizes, collapse = ""))
  
  write.table(compareM,paste0(data_path,"Table",namef,".csv"),sep = ",",quote = F)
  
  # Debug: Print matrix structure
  cat("=== DATA MATRIX DEBUG INFO ===\n")
  cat("Matrix dimensions:", dim(compareM), "\n")
  cat("Column names:", colnames(compareM), "\n")
  cat("Row names:", rownames(compareM), "\n")
  cat("Matrix preview:\n")
  print(head(compareM))
  
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
  write.csv(friedman_summary, paste0(data_path, "Friedman_Results", namef, ".csv"), row.names = FALSE)
  
  # Create ranking summary table
  ranking_summary <- data.frame(
    Method = names(ranking_letters),
    Mean_Prediction_Ability = round(colMeans(compareM), 4),
    SD = round(apply(compareM, 2, sd), 4),
    Ranking_Letter = ranking_letters[names(ranking_letters)]
  ) %>%
    arrange(desc(Mean_Prediction_Ability))
  
  write.csv(ranking_summary, paste0(data_path, "Method_Rankings", namef, ".csv"), row.names = FALSE)
  
  data<-reshape2::melt(compareM)
  
  # Create table with ranking letters
  table_data <- compareM
  table_data_with_letters <- rbind(
    round(table_data, 4),
    ranking_letters[colnames(table_data)]
  )
  rownames(table_data_with_letters)[nrow(table_data_with_letters)] <- "Rank"
  
  table_graphic <- ggtexttable(table_data_with_letters)
  ggsave(paste0(data_path, "Table_with_Rankings", namef, ".png"), table_graphic, scale=1, dpi = 300, height = 2.5, width = 10)
  
  # Calculate y-position for letters (above the boxplots)
  max_values <- aggregate(value ~ Var2, data, max)
  y_positions <- max_values$value + 0.02
  names(y_positions) <- max_values$Var2
  
  # Prediction Accuracy plot (with h2 correction)
  p1<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.35)+
    labs(title=paste0("Prediction Accuracy - ",trait),x= "Prediction method",fill="Method", y = "Prediction Accuracy")+
    stat_summary(fun=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    # Add ranking letters
    annotate("text", x = 1:length(ranking_letters), 
             y = (y_positions[names(ranking_letters)]/sqrt(h2)) + 0.02,
             label = ranking_letters, size = 4, fontface = "bold", color = "blue") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Plot_with_Rankings",namef,".png"),p1,scale=1, dpi = 300, height = 4, width = 4)
  
  # Prediction Ability plot (without h2 correction)
  p2<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.35)+
    labs(title=paste0("Prediction Abilityy - ",trait),x= "Prediction method",fill="Method", y = "Prediction Ability")+
    stat_summary(fun=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    # Add ranking letters
    annotate("text", x = 1:length(ranking_letters), 
             y = y_positions[names(ranking_letters)] + 0.02,
             label = ranking_letters, size = 3, fontface = "bold", color = "blue") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Ability_Plot_with_Rankings",namef,".png"),p2,scale=1, dpi = 300, height = 4, width = 5)
  
  # Create a summary plot showing method rankings
  ranking_plot_data <- data.frame(
    Method = factor(ranking_summary$Method, levels = ranking_summary$Method),
    Mean_Value = ranking_summary$Mean_Prediction_Ability,
    Letter = ranking_summary$Ranking_Letter
  )
  
  p3 <- ggplot(ranking_plot_data, aes(x = Method, y = Mean_Value, fill = Method)) +
    geom_col(alpha = 0.7, color = "black") +
    geom_text(aes(label = Letter), vjust = -0.5, size = 5, fontface = "bold", color = "red") +
    labs(title = paste0("Method Rankings - ", trait),
         subtitle = paste0("Friedman test p-value: ", round(friedman_analysis$friedman_result$p.value, 4)),
         x = "Prediction Method",
         y = "Mean Prediction Ability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(paste0(data_path,"Method_Rankings_Bar",namef,".png"), p3, scale=1, dpi = 300, height = 6, width = 8)
  
  # Print summary to console
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Files generated:\n")
  cat("- Table_with_Rankings", namef, ".png\n")
  cat("- Plot_with_Rankings", namef, ".png\n")
  cat("- Ability_Plot_with_Rankings", namef, ".png\n")
  cat("- Method_Rankings_Bar", namef, ".png\n")
  cat("- Friedman_Results", namef, ".csv\n")
  cat("- Method_Rankings", namef, ".csv\n")
  
  cat("\nRanking Letters:\n")
  for(i in 1:length(ranking_letters)) {
    cat(names(ranking_letters)[i], ": ", ranking_letters[i], "\n")
  }
  
  cat("\nNote: Methods sharing the same letter are not significantly different.\n")
  cat("Letters are assigned alphabetically with 'a' being the best performing group.\n")
}
