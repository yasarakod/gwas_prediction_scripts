# ============================================================================== #
# Script: 12_multivariate_genetic_correlation.R
# Title: Multivariate Genetic & Phenotypic Correlation Analysis
# Description:
#   Fits multivariate mixed models to calculate genetic, phenotypic, and environmental correlation matrices among survival and morphological traits, outputting correlation plots.
#
# Inputs:
#   - data/pheno (Multi-trait phenotype records)
#   - data/M (Genotype matrix)
#   - data/A_matrix.rds (Pedigree matrix)
#
# Outputs:
#   - results/genetic_correlations_all_methods.rds (Correlation matrices)
#   - results/genetic_correlations_gblup.png (Visual correlation correlogram)
#
# Required R Packages: rrBLUP, BGLR, sommer, AGHmatrix, corrplot, ggplot2
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


# Genetic Correlation Analysis for Multiple Traits
# setwd() configured via project environment

# Load required libraries
library(rrBLUP)
library(BGLR)
library(sommer)  # For multivariate models
library(AGHmatrix)
library(corrplot)
library(ggplot2)
library(reshape2)

# Method 1: Multivariate GBLUP for Genetic Correlations
calculate_genetic_correlations_gblup <- function(trait_names = c("trait1", "trait2", "trait3")) {
  
  cat("=== Multivariate GBLUP Genetic Correlation Analysis ===\n")
  
  # Load data
  M <- readRDS("data/M")
  pheno <- readRDS("data/pheno")
  
  cat("Available traits:", colnames(pheno), "\n")
  cat("Analyzing traits:", paste(trait_names, collapse = ", "), "\n")
  
  # Check if all traits exist
  missing_traits <- trait_names[!trait_names %in% colnames(pheno)]
  if(length(missing_traits) > 0) {
    stop("Traits not found: ", paste(missing_traits, collapse = ", "))
  }
  
  # Extract trait data
  Y <- pheno[, trait_names, drop = FALSE]
  
  # Remove individuals with missing data for any trait
  complete_cases <- complete.cases(Y)
  Y_clean <- Y[complete_cases, ]
  M_clean <- M[complete_cases, ]
  
  cat("Complete observations across all traits:", nrow(Y_clean), "\n")
  cat("Number of markers:", ncol(M_clean), "\n")
  
  # Create genomic relationship matrix
  cat("Creating genomic relationship matrix...\n")
  G <- A.mat(M_clean - 1)
  
  # Fit multivariate GBLUP using sommer
  cat("Fitting multivariate GBLUP model...\n")
  
  tryCatch({
    # Convert to matrix format for sommer
    Y_matrix <- as.matrix(Y_clean)
    
    # Create data frame for sommer
    data_mv <- data.frame(
      id = rownames(Y_clean),
      Y_matrix
    )
    
    # Fit multivariate model
    model_mv <- mmer(
      fixed = cbind(get(trait_names[1]), get(trait_names[2]), get(trait_names[3])) ~ 1,
      random = ~ vs(id, Gu = G),
      rcov = ~ vs(units),
      data = data_mv,
      verbose = FALSE
    )
    
    # Extract genetic variance-covariance matrix
    G_var <- model_mv$sigma$`u:id`
    
    # Extract residual variance-covariance matrix  
    R_var <- model_mv$sigma$units
    
    # Calculate genetic correlations
    G_corr <- cov2cor(G_var)
    
    # Calculate phenotypic correlations
    P_corr <- cor(Y_clean)
    
    # Print results
    cat("\n=== GENETIC VARIANCE-COVARIANCE MATRIX ===\n")
    print(round(G_var, 6))
    
    cat("\n=== GENETIC CORRELATION MATRIX ===\n")
    print(round(G_corr, 3))
    
    cat("\n=== PHENOTYPIC CORRELATION MATRIX ===\n")
    print(round(P_corr, 3))
    
    # Calculate heritabilities
    P_var <- G_var + R_var
    h2_diag <- diag(G_var) / diag(P_var)
    names(h2_diag) <- trait_names
    
    cat("\n=== HERITABILITIES ===\n")
    for(i in 1:length(h2_diag)) {
      cat(names(h2_diag)[i], ": h² =", round(h2_diag[i], 3), "\n")
    }
    
    # Create results object
    results <- list(
      method = "Multivariate GBLUP",
      traits = trait_names,
      n_obs = nrow(Y_clean),
      n_markers = ncol(M_clean),
      genetic_correlations = G_corr,
      phenotypic_correlations = P_corr,
      genetic_variances = diag(G_var),
      heritabilities = h2_diag,
      genetic_var_cov = G_var,
      residual_var_cov = R_var,
      model = model_mv
    )
    
    return(results)
    
  }, error = function(e) {
    cat("Error in multivariate GBLUP:", e$message, "\n")
    cat("Trying alternative approach with BGLR...\n")
    return(calculate_genetic_correlations_bglr(trait_names, Y_clean, G))
  })
}

# Alternative method using BGLR for genetic correlations
calculate_genetic_correlations_bglr <- function(trait_names, Y_clean, G) {
  
  cat("=== Alternative: BGLR-based Genetic Correlations ===\n")
  
  # Fit univariate models for each trait
  models <- list()
  breeding_values <- matrix(NA, nrow = nrow(Y_clean), ncol = length(trait_names))
  colnames(breeding_values) <- trait_names
  rownames(breeding_values) <- rownames(Y_clean)
  
  for(i in 1:length(trait_names)) {
    cat("Fitting model for trait:", trait_names[i], "\n")
    
    y <- Y_clean[, i]
    
    # Set up BGLR
    ETA <- list(list(K = G, model = "RKHS"))
    
    # Fit model
    model <- BGLR(y = y, ETA = ETA, nIter = 5000, burnIn = 1000, verbose = FALSE)
    
    models[[trait_names[i]]] <- model
    breeding_values[, i] <- model$ETA[[1]]$u
  }
  
  # Calculate genetic correlations from breeding values
  genetic_correlations <- cor(breeding_values)
  phenotypic_correlations <- cor(Y_clean)
  
  cat("\n=== GENETIC CORRELATIONS (from breeding values) ===\n")
  print(round(genetic_correlations, 3))
  
  results <- list(
    method = "BGLR univariate",
    traits = trait_names,
    genetic_correlations = genetic_correlations,
    phenotypic_correlations = phenotypic_correlations,
    breeding_values = breeding_values,
    models = models
  )
  
  return(results)
}

# Method 2: PBLUP-based genetic correlations
calculate_genetic_correlations_pblup <- function(trait_names = c("trait1", "trait2", "trait3")) {
  
  cat("=== PBLUP Genetic Correlation Analysis ===\n")
  
  # Load data
  pheno <- readRDS("data/pheno")
  
  # Check if A matrix exists
  if(!file.exists("data/A_matrix.rds")) {
    stop("A matrix not found.Please run pedigree processing first.")
  }
  
  A_matrix <- readRDS("data/A_matrix.rds")
  
  # Match individuals
  common_ids <- intersect(pheno$Sample_ID, rownames(A_matrix))
  cat("Common individuals:", length(common_ids), "\n")
  
  # Extract trait data for common individuals
  pheno_matched <- pheno[pheno$Sample_ID %in% common_ids, ]
  A_sub <- A_matrix[common_ids, common_ids]
  Y <- pheno_matched[, trait_names, drop = FALSE]
  
  # Remove individuals with missing data
  complete_cases <- complete.cases(Y)
  Y_clean <- Y[complete_cases, ]
  A_clean <- A_sub[complete_cases, complete_cases]
  
  cat("Complete observations:", nrow(Y_clean), "\n")
  
  # Fit multivariate PBLUP
  tryCatch({
    # Create data for sommer
    data_mv <- data.frame(
      id = rownames(Y_clean),
      Y_clean
    )
    
    # Fit multivariate model
    model_mv <- mmer(
      fixed = cbind(get(trait_names[1]), get(trait_names[2]), get(trait_names[3])) ~ 1,
      random = ~ vs(id, Gu = A_clean),
      rcov = ~ vs(units),
      data = data_mv,
      verbose = FALSE
    )
    
    # Extract results (same as GBLUP)
    G_var <- model_mv$sigma$`u:id`
    R_var <- model_mv$sigma$units
    G_corr <- cov2cor(G_var)
    P_corr <- cor(Y_clean)
    
    cat("\n=== PBLUP GENETIC CORRELATIONS ===\n")
    print(round(G_corr, 3))
    
    results <- list(
      method = "Multivariate PBLUP",
      traits = trait_names,
      genetic_correlations = G_corr,
      phenotypic_correlations = P_corr,
      genetic_var_cov = G_var,
      residual_var_cov = R_var
    )
    
    return(results)
    
  }, error = function(e) {
    cat("Error in multivariate PBLUP:", e$message, "\n")
    return(NULL)
  })
}

# Method 3: Compare genetic correlations between methods
compare_genetic_correlations <- function(trait_names = c("trait1", "trait2", "trait3")) {
  
  cat("=== Comparing Genetic Correlations: GBLUP vs PBLUP ===\n")
  
  # Calculate GBLUP correlations
  gblup_results <- calculate_genetic_correlations_gblup(trait_names)
  
  # Calculate PBLUP correlations
  pblup_results <- calculate_genetic_correlations_pblup(trait_names)
  
  if(!is.null(gblup_results) && !is.null(pblup_results)) {
    
    # Extract correlation matrices
    G_corr_gblup <- gblup_results$genetic_correlations
    G_corr_pblup <- pblup_results$genetic_correlations
    
    # Print comparison
    cat("\n=== COMPARISON SUMMARY ===\n")
    cat("GBLUP Genetic Correlations:\n")
    print(round(G_corr_gblup, 3))
    
    cat("\nPBLUP Genetic Correlations:\n")
    print(round(G_corr_pblup, 3))
    
    cat("\nDifferences (GBLUP - PBLUP):\n")
    print(round(G_corr_gblup - G_corr_pblup, 3))
    
    # Create comparison object
    comparison <- list(
      traits = trait_names,
      gblup_correlations = G_corr_gblup,
      pblup_correlations = G_corr_pblup,
      differences = G_corr_gblup - G_corr_pblup,
      gblup_results = gblup_results,
      pblup_results = pblup_results
    )
    
    return(comparison)
  }
  
  return(list(gblup_results = gblup_results, pblup_results = pblup_results))
}

# Visualization function
plot_genetic_correlations <- function(results, title = "Genetic Correlations") {
  
  if("genetic_correlations" %in% names(results)) {
    corr_matrix <- results$genetic_correlations
  } else {
    stop("Results object must contain genetic_correlations")
  }
  
  # Create correlation plot
  corrplot(corr_matrix, 
           method = "color",
           type = "upper",
           order = "original",
           tl.cex = 1.2,
           tl.col = "black",
           addCoef.col = "black",
           number.cex = 1.2,
           title = title,
           mar = c(0,0,2,0))
  
  # Also create ggplot version
  corr_melted <- melt(corr_matrix)
  
  p <- ggplot(corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Genetic\nCorrelation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = title, x = "", y = "") +
    coord_fixed()
  
  return(p)
}

# Usage examples and data check
check_genetic_correlation_data <- function() {
  cat("=== GENETIC CORRELATION DATA CHECK ===\n")
  
  files_needed <- c("data/M.rds", "data/pheno.rds", "data/A_matrix.rds")
  
  for(file in files_needed) {
    if(file.exists(file)) {
      cat("✓", file, "- Available\n")
    } else {
      cat("✗", file, "- Missing\n")
    }
  }
  
  # Check phenotype data
  if(file.exists("data/pheno.rds")) {
    pheno <- readRDS("data/pheno.rds")
    cat("\nAvailable traits:", paste(colnames(pheno), collapse = ", "), "\n")
    
    # Check for your specific traits
    your_traits <- c("Is_Dead", "DPC_date", "DPC_time")
    available_traits <- your_traits[your_traits %in% colnames(pheno)]
    missing_traits <- your_traits[!your_traits %in% colnames(pheno)]
    
    if(length(available_traits) > 0) {
      cat("Your traits available:", paste(available_traits, collapse = ", "), "\n")
    }
    if(length(missing_traits) > 0) {
      cat("Your traits missing:", paste(missing_traits, collapse = ", "), "\n")
    }
  }
}

# Main analysis function
run_genetic_correlation_analysis <- function() {
  
  # Your specific traits
  traits <- c("Is_Dead", "DPC_date", "DPC_time")
  
  cat("=== GENETIC CORRELATION ANALYSIS FOR YOUR TRAITS ===\n")
  cat("Traits:", paste(traits, collapse = ", "), "\n\n")
  
  # Run GBLUP analysis
  cat("1.GBLUP Analysis:\n")
  gblup_results <- calculate_genetic_correlations_gblup(traits)
  
  # Run PBLUP analysis (if A matrix available)
  if(file.exists("data/A_matrix.rds")) {
    cat("\n2.PBLUP Analysis:\n")
    pblup_results <- calculate_genetic_correlations_pblup(traits)
    
    # Compare methods
    cat("\n3.Method Comparison:\n")
    comparison <- compare_genetic_correlations(traits)
    
    # Save all results
    saveRDS(list(gblup = gblup_results, pblup = pblup_results, comparison = comparison),
            "data/genetic_correlations_all_methods.rds")
    
  } else {
    cat("\nPBLUP analysis skipped - A matrix not available\n")
    saveRDS(gblup_results, "data/genetic_correlations_gblup.rds")
  }
  
  # Create visualizations
  if(exists("gblup_results") && !is.null(gblup_results)) {
    png("genetic_correlations_gblup.png", width = 800, height = 600)
    plot_genetic_correlations(gblup_results, "GBLUP Genetic Correlations")
    dev.off()
    cat("GBLUP correlation plot saved as genetic_correlations_gblup.png\n")
  }
  
  return(gblup_results)
}

# Display usage instructions
cat("=== USAGE INSTRUCTIONS ===\n")
cat("1.Check data availability:\n")
cat("   check_genetic_correlation_data()\n\n")

cat("2.Run full analysis for your traits:\n")
cat("   results <- run_genetic_correlation_analysis()\n\n")

cat("3.Custom trait analysis:\n")
cat("   results <- calculate_genetic_correlations_gblup(c('trait1', 'trait2', 'trait3'))\n\n")

cat("4.Compare methods:\n")
cat("   comparison <- compare_genetic_correlations(c('trait1', 'trait2', 'trait3'))\n\n")

# Run data check
check_genetic_correlation_data()

# Run complete analysis for your three traits
results <- run_genetic_correlation_analysis()
