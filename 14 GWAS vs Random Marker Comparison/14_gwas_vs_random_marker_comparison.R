# ============================================================================== #
# Script: 14_gwas_vs_random_marker_comparison.R
# Title: GWAS vs Random Marker Statistical Comparison & Plotting
# Description:
#   Analyzes cross-validation prediction results across varying marker densities (50 to ALL), performing Friedman non-parametric tests and Nemenyi post-hoc ranking letter assignments.
#
# Inputs:
#   - results/GBLUP_Folds_5_Times_10_* (Cross-validation output text files for random and top markers)
#
# Outputs:
#   - results/Table_<pmethod><namef>.csv (Statistical summary table)
#   - results/Friedman_Results_<pmethod><namef>.csv (Friedman test statistics and p-values)
#   - results/Marker_Rankings_<pmethod><namef>.csv (Rankings and significance grouping letters)
#   - results/Accuracy_Plot_with_Rankings_*.png (Accuracy boxplot with ranking letters)
#   - results/Ability_Plot_with_Rankings_*.png (Predictive ability boxplot with ranking letters)
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


#Plotting GWAS and Random Marker Comparison Results with Friedman Test and Ranking Letters
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

# Parse Arguments
dir<-"data"
graph_type<-"B"
trait_name<-c("Is_Dead")
pmethods<-c("GBLUP")
iter<-10
fold<-5
pop_sizes<-"NULL"
Fixed_effect<-c("NULL")
h2<-1

# Define two marker sets
random_markers <- c(50, 100, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000, "NULL")
gwas_markers <- c("top_50", "top_100", "top_500", "top_1000", "top_5000", "top_10000", "top_20000", "top_30000", "top_40000", "top_50000")

# Combine all markers for processing
all_markers <- c(random_markers, gwas_markers)

data_path = paste0(dir, "/")

if (graph_type=="B"){
  
  pmethod<-pmethods
  sname<-trait_name
  
  # Create comparison matrix
  compareM <- matrix(0, iter, length(all_markers))
  
  for (i in seq_along(all_markers)) {  
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, "_", sname, "_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_", all_markers[i], ".txt")
    
    # Check if file exists
    if (file.exists(filename)) {
      datap <- read.table(filename)
      for (n in 1:iter) {
        compareM[n, i] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
      }
    } else {
      cat("Warning: File not found:", filename, "\n")
      compareM[, i] <- NA  # Fill with NA if file doesn't exist
    }
  }
  
  # Set column names with 'm=' prefix
  colnames(compareM) <- paste0("m=", all_markers)
  rownames(compareM) <- paste0("CV", 1:iter)
  
  # Remove columns with all NA values
  compareM <- compareM[, !apply(is.na(compareM), 2, all)]
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_",sname,"_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_marker_comparison")
  
  # Save the comparison matrix
  write.table(compareM, paste0(data_path,"Table_",pmethod,namef,".csv"), sep = ",", quote = FALSE)
  
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
    Marker_Config = names(ranking_letters),
    Mean_Prediction_Ability = round(colMeans(compareM, na.rm = TRUE), 4),
    SD = round(apply(compareM, 2, sd, na.rm = TRUE), 4),
    Ranking_Letter = ranking_letters[names(ranking_letters)]
  ) %>%
    arrange(desc(Mean_Prediction_Ability))
  
  write.csv(ranking_summary, paste0(data_path, "Marker_Rankings_", pmethod, namef, ".csv"), row.names = FALSE)
  
  # Melt the data for plotting
  plot_data <- reshape2::melt(compareM)
  
  # Extract numeric values and create grouping
  plot_data$numeric_value <- as.numeric(gsub(".*?([0-9]+).*", "\\1", plot_data$Var2))
  plot_data$numeric_value[grepl("NULL", plot_data$Var2)] <- "NULL"  # Handle NULL case
  plot_data$group <- ifelse(grepl("top_", plot_data$Var2), "GWAS", "Random")
  
  # Create proper ordering for x-axis
  marker_order <- c("50", "100", "500", "1000", "5000", "10000", "20000", "30000", "40000", "50000", "NULL")
  plot_data$x_label <- factor(plot_data$numeric_value, levels = marker_order)
  
  # Calculate y-positions for ranking letters
  max_values <- plot_data %>%
    group_by(x_label, group) %>%
    summarise(max_val = max(value, na.rm = TRUE), .groups = 'drop')
  
  # Add ranking letters to max_values
  max_values$ranking_letter <- ""
  for (i in 1:nrow(max_values)) {
    # Find corresponding marker configuration
    marker_val <- as.character(max_values$x_label[i])
    group_val <- max_values$group[i]
    
    if (group_val == "GWAS" && marker_val != "NULL") {
      config_name <- paste0("m=top_", marker_val)
    } else if (group_val == "Random") {
      config_name <- paste0("m=", marker_val)
    }
    
    if (exists("config_name") && config_name %in% names(ranking_letters)) {
      max_values$ranking_letter[i] <- ranking_letters[config_name]
    }
  }
  
  # Create the first plot - Prediction Accuracy (with h2 correction)
  p1 <- ggplot(plot_data, aes(x=x_label, y=value/sqrt(h2), fill=group)) + 
    geom_boxplot(position = position_dodge(width = 0.8), color="black", alpha = 0.8) +
    scale_fill_manual(
      values = c("Random" = "#443A83FF", "GWAS" = "#27AD81FF"),
      name = "Marker Type"
    ) +
    labs(title=paste0("Prediction Accuracy - ", sname), 
         x="Number of Markers", 
         y="Prediction Accuracy") +
    stat_summary(aes(group = group), fun=mean, geom="point", shape=21, size=2, 
                 color="black", fill="red", stroke=1,
                 position = position_dodge(width = 0.8)) +
    # Add ranking letters
    geom_text(data = max_values, 
              aes(x = x_label, y = (max_val/sqrt(h2)) + 0.02, 
                  label = ranking_letter, group = group),
              position = position_dodge(width = 0.8),
              size = 4, fontface = "bold", color = "blue", inherit.aes = FALSE) +
    theme(legend.title = element_text(colour="blue", size=12, face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Accuracy_Plot_with_Rankings_",pmethod,namef,".png"), p1, scale=1, dpi = 300, height = 6, width = 12)
  
  # Create the second plot - Predictive Ability (without h2 correction)
  p2 <- ggplot(plot_data, aes(x=x_label, y=value, fill=group)) + 
    geom_boxplot(position = position_dodge(width = 0.8), color="black", alpha = 0.8) +
    scale_fill_manual(
      values = c("Random" = "#443A83FF", "GWAS" = "#27AD81FF"),
      name = "Marker Type"
    ) +
    labs(title=paste0(sname), 
         x="Number of Markers", 
         y="Predictive Ability") +
    stat_summary(aes(group = group), fun=mean, geom="point", shape=21, size=2, 
                 color="black", fill="red", stroke=1,
                 position = position_dodge(width = 0.8)) +
    # Add ranking letters
    geom_text(data = max_values, 
              aes(x = x_label, y = max_val + 0.02, 
                  label = ranking_letter, group = group),
              position = position_dodge(width = 0.8),
              size = 5, fontface = "bold", color = "blue", inherit.aes = FALSE) +
    theme(
      legend.title = element_text(colour="blue", size=16, face="bold", margin=margin(b=15)),  # Increased legend title size and added margin
      legend.text = element_text(size=14, face="bold"),  # Added bold legend text
      plot.title = element_text(hjust=0.5, size=18, face="bold"),  # Made title bold
      axis.title = element_text(size=18, face="bold", margin=margin(t=20)),  # Made axis titles bold
      axis.text = element_text(size=15, face="bold"),   # Made axis text bold
      axis.text.x = element_text(angle=45, hjust=1, size=15, face="bold"),  # Bold x-axis text
      axis.text.y = element_text(size=15, face="bold"),  # Bold y-axis text
      plot.background = element_rect(fill = "#FFFFFF"),
      legend.background = element_rect(linetype = "dashed"),
      panel.border = element_rect(fill = NA)
    )
  
  ggsave(paste0(data_path,"Ability_Plot_with_Rankings_",pmethod,namef,".png"), 
         p2, scale=1, dpi = 300, height = 6, width = 12)
  # Print summary to console
  cat("\n=== MARKER COMPARISON ANALYSIS COMPLETE ===\n")
  cat("Files generated:\n")
  cat("- Accuracy_Plot_with_Rankings_", pmethod, namef, ".png\n")
  cat("- Ability_Plot_with_Rankings_", pmethod, namef, ".png\n")
  
  cat("\nMarker Configuration Ranking Letters:\n")
  for(i in 1:length(ranking_letters)) {
    cat(names(ranking_letters)[i], ": ", ranking_letters[i], "\n")
  }
  
  cat("\nNote: Configurations sharing the same letter are not significantly different.\n")
  cat("Letters are assigned alphabetically with 'a' being the best performing group.\n")
}
