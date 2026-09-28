# ============================================================================== #
# Script: 11_pblup_vs_gblup_heritability.R
# Title: PBLUP vs GBLUP Heritability Comparison
# Description:
#   Estimates pedigree-based heritability (PBLUP) using the additive relationship matrix (A) and performs comparative evaluation against genomic heritability (GBLUP).
#
# Inputs:
#   - data/pheno (Phenotype observations)
#   - data/A_matrix.rds (Numerator relationship matrix)
#   - data/M (Genotype matrix for GBLUP comparison)
#
# Outputs:
#   - results/pblup_heritability_<trait>.rds (PBLUP heritability results)
#   - results/gblup_pblup_comparison_<trait>.rds (Model comparison statistics)
#   - results/gblup_pblup_heritability_summary.csv (Tabular summary of h2 estimates)
#
# Required R Packages: rrBLUP, AGHmatrix
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


# PBLUP Heritability Calculator from Your Processed Data Files
# Set working directory to your data path
# setwd() configured via project environment

# Load required libraries
library(rrBLUP)
library(AGHmatrix)  # For pedigree-based relationship matrix

# Method 1: PBLUP using pre-computed A matrix (if available)
calculate_pblup_heritability_with_Amatrix <- function(trait_column = "trait_name") {
  
  cat("=== PBLUP Heritability Analysis ===\n")
  
  # Load processed phenotype data
  pheno <- readRDS("data/pheno")
  cat("Phenotype data dimensions:", dim(pheno), "\n")
  cat("Available traits:", colnames(pheno), "\n")
  
  # Check if A matrix exists from previous analysis
  if(file.exists("data/A_matrix.rds")) {
    cat("Loading pre-computed A matrix...\n")
    A_matrix <- readRDS("data/A_matrix.rds")
    cat("A matrix dimensions:", dim(A_matrix), "\n")
  } else {
    stop("A matrix not found.Please run pedigree processing first or use Method 2.")
  }
  
  # Select trait for analysis
  if(trait_column %in% colnames(pheno)) {
    y <- pheno[, trait_column]
    cat("Analyzing trait:", trait_column, "\n")
  } else {
    cat("Available traits:", paste(colnames(pheno), collapse = ", "), "\n")
    stop("Trait '", trait_column, "' not found in phenotype data")
  }
  
  # Match individuals between phenotype and A matrix
  common_ids <- intersect(pheno$Sample_ID, rownames(A_matrix))
  cat("Common individuals between phenotype and A matrix:", length(common_ids), "\n")
  
  if(length(common_ids) == 0) {
    stop("No common individuals found between phenotype and A matrix!")
  }
  
  # Subset data to common individuals
  pheno_matched <- pheno[pheno$Sample_ID %in% common_ids, ]
  A_sub <- A_matrix[common_ids, common_ids]
  y_matched <- pheno_matched[, trait_column]
  names(y_matched) <- pheno_matched$Sample_ID
  
  # Remove individuals with missing phenotypes
  complete_cases <- !is.na(y_matched)
  y_clean <- y_matched[complete_cases]
  A_clean <- A_sub[complete_cases, complete_cases]
  
  cat("Complete observations for PBLUP:", length(y_clean), "\n")
  cat("Trait range:", round(range(y_clean), 3), "\n")
  
  # Fit PBLUP model
  cat("Fitting PBLUP model...\n")
  model <- mixed.solve(y = y_clean, K = A_clean)
  
  # Calculate heritability
  Vg <- model$Vu  # Additive genetic variance
  Ve <- model$Ve  # Error variance
  Vp <- Vg + Ve   # Phenotypic variance
  h2 <- Vg / Vp   # Narrow-sense heritability
  
  # Standard error
  n_obs <- length(y_clean)
  se_h2 <- sqrt((4 * h2 * (1 - h2)^2) / n_obs)
  
  # Print results
  cat("\n=== PBLUP HERITABILITY RESULTS ===\n")
  cat("Method: Pedigree-Based BLUP\n")
  cat("Trait:", trait_column, "\n")
  cat("Sample size:", n_obs, "\n")
  cat("h² =", round(h2, 3), "± SE =", round(se_h2, 3), "\n")
  cat("Additive genetic variance (Va) =", sprintf("%.6f", Vg), "\n")
  cat("Error variance (Ve) =", sprintf("%.6f", Ve), "\n")
  cat("Phenotypic variance (Vp) =", sprintf("%.6f", Vp), "\n")
  cat("Mean A matrix diagonal =", round(mean(diag(A_clean)), 3), "\n")
  cat("Mean A matrix off-diagonal =", round(mean(A_clean[upper.tri(A_clean)]), 3), "\n")
  
  # Save results
  results <- list(
    method = "PBLUP",
    trait = trait_column,
    h2 = h2,
    se_h2 = se_h2,
    Va = Vg,  # Additive genetic variance
    Ve = Ve,
    Vp = Vp,
    n_obs = n_obs,
    A_matrix_stats = list(
      diagonal_mean = mean(diag(A_clean)),
      offdiagonal_mean = mean(A_clean[upper.tri(A_clean)])
    )
  )
  
  saveRDS(results, paste0("data/pblup_heritability_", trait_column, ".rds"))
  cat("Results saved to: data/pblup_heritability_", trait_column, ".rds\n")
  
  return(results)
}

# Method 2: PBLUP using pedigree file (if A matrix not available)
calculate_pblup_heritability_from_pedigree <- function(pedigree_file = "data/pedigree.txt", 
                                                       trait_column = "trait_name") {
  
  cat("=== PBLUP Heritability from Pedigree File ===\n")
  
  # Check if pedigree file exists
  if(!file.exists(pedigree_file)) {
    stop("Pedigree file not found at ", pedigree_file)
  }
  
  # Load pedigree data
  cat("Reading pedigree data from", pedigree_file, "\n")
  ped_df <- read.table(pedigree_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
  
  # Ensure column names are correct
  colnames(ped_df)[1:3] <- c("ID", "Sire", "Dam")
  cat("Original pedigree has", nrow(ped_df), "animals\n")
  
  # Clean pedigree data
  ped_df$Sire[is.na(ped_df$Sire) | ped_df$Sire == ""] <- "0"
  ped_df$Dam[is.na(ped_df$Dam) | ped_df$Dam == ""] <- "0"
  
  # Add missing parents to pedigree
  all_ids <- unique(ped_df$ID)
  all_sires <- unique(ped_df$Sire)
  all_dams <- unique(ped_df$Dam)
  
  missing_sires <- setdiff(all_sires[all_sires != "0"], all_ids)
  missing_dams <- setdiff(all_dams[all_dams != "0"], all_ids)
  all_missing <- unique(c(missing_sires, missing_dams))
  
  if(length(all_missing) > 0) {
    missing_rows <- data.frame(
      ID = all_missing,
      Sire = "0",
      Dam = "0",
      stringsAsFactors = FALSE
    )
    ped_df <- rbind(ped_df, missing_rows)
  }
  
  cat("Completed pedigree has", nrow(ped_df), "animals\n")
  
  # Create A matrix
  cat("Creating additive relationship matrix (A)...\n")
  A_matrix <- Amatrix(ped_df[,1:3])
  
  # Save A matrix for future use
  saveRDS(A_matrix, "data/A_matrix.rds")
  saveRDS(ped_df, "data/completed_pedigree.rds")
  cat("A matrix saved for future use\n")
  
  # Now calculate heritability using the A matrix
  calculate_pblup_heritability_with_Amatrix(trait_column)
}

# Method 3: Compare GBLUP vs PBLUP heritability
compare_gblup_pblup_heritability <- function(trait_column = "trait_name") {
  
  cat("=== GBLUP vs PBLUP Heritability Comparison ===\n")
  
  # Calculate GBLUP heritability
  cat("Calculating GBLUP heritability...\n")
  M <- readRDS("data/M")
  pheno <- readRDS("data/pheno")
  
  # GBLUP analysis
  y_gblup <- pheno[, trait_column]
  complete_cases_g <- !is.na(y_gblup)
  M_clean <- M[complete_cases_g, ]
  y_gblup_clean <- y_gblup[complete_cases_g]
  
  G <- A.mat(M_clean - 1)
  model_gblup <- mixed.solve(y = y_gblup_clean, K = G)
  
  h2_gblup <- model_gblup$Vu / (model_gblup$Vu + model_gblup$Ve)
  
  # Calculate PBLUP heritability
  cat("Calculating PBLUP heritability...\n")
  pblup_results <- calculate_pblup_heritability_with_Amatrix(trait_column)
  h2_pblup <- pblup_results$h2
  
  # Comparison results
  cat("\n=== COMPARISON RESULTS ===\n")
  cat("Trait:", trait_column, "\n")
  cat("GBLUP h² =", round(h2_gblup, 3), "\n")
  cat("PBLUP h² =", round(h2_pblup, 3), "\n")
  cat("Difference (GBLUP - PBLUP) =", round(h2_gblup - h2_pblup, 3), "\n")
  
  if(h2_gblup > h2_pblup) {
    cat("GBLUP captures", round((h2_gblup - h2_pblup) / h2_pblup * 100, 1), "% more genetic variance than PBLUP\n")
  } else if(h2_pblup > h2_gblup) {
    cat("PBLUP captures", round((h2_pblup - h2_gblup) / h2_gblup * 100, 1), "% more genetic variance than GBLUP\n")
  } else {
    cat("Both methods estimate similar heritability\n")
  }
  
  comparison <- list(
    trait = trait_column,
    h2_gblup = h2_gblup,
    h2_pblup = h2_pblup,
    difference = h2_gblup - h2_pblup,
    gblup_results = list(Vg = model_gblup$Vu, Ve = model_gblup$Ve),
    pblup_results = pblup_results
  )
  
  saveRDS(comparison, paste0("data/gblup_pblup_comparison_", trait_column, ".rds"))
  return(comparison)
}

# Method 4: Analyze all traits with both methods
analyze_all_traits_both_methods <- function() {
  
  pheno <- readRDS("data/pheno")
  trait_columns <- names(pheno)[sapply(pheno, is.numeric)]
  trait_columns <- trait_columns[trait_columns != "Sample_ID"]
  
  cat("Analyzing", length(trait_columns), "traits with both GBLUP and PBLUP\n\n")
  
  results_summary <- data.frame(
    Trait = character(),
    GBLUP_h2 = numeric(),
    PBLUP_h2 = numeric(),
    Difference = numeric(),
    GBLUP_n = integer(),
    PBLUP_n = integer(),
    stringsAsFactors = FALSE
  )
  
  for(trait in trait_columns) {
    cat("=== Analyzing", trait, "===\n")
    tryCatch({
      comparison <- compare_gblup_pblup_heritability(trait)
      
      results_summary <- rbind(results_summary, data.frame(
        Trait = trait,
        GBLUP_h2 = round(comparison$h2_gblup, 3),
        PBLUP_h2 = round(comparison$h2_pblup, 3),
        Difference = round(comparison$difference, 3),
        GBLUP_n = length(comparison$gblup_results$Vg),
        PBLUP_n = comparison$pblup_results$n_obs
      ))
      
    }, error = function(e) {
      cat("Error analyzing", trait, ":", e$message, "\n")
    })
    cat("\n")
  }
  
  cat("=== SUMMARY: GBLUP vs PBLUP ===\n")
  print(results_summary)
  
  write.csv(results_summary, "data/gblup_pblup_heritability_summary.csv", row.names = FALSE)
  return(results_summary)
}

# Check data availability for PBLUP
check_pblup_data_availability <- function() {
  cat("=== PBLUP DATA AVAILABILITY CHECK ===\n")
  
  files_to_check <- c(
    "data/pheno.rds",
    "data/A_matrix.rds", 
    "data/pedigree.txt",
    "data/completed_pedigree.rds"
  )
  
  for(file in files_to_check) {
    if(file.exists(file)) {
      cat("✓", file, "- Available\n")
    } else {
      cat("✗", file, "- Missing\n")
    }
  }
  
  cat("\nRecommended usage order:\n")
  cat("1.If A_matrix.rds exists: calculate_pblup_heritability_with_Amatrix('trait_name')\n")
  cat("2.If only pedigree.txt exists: calculate_pblup_heritability_from_pedigree('data/pedigree.txt', 'trait_name')\n")
  cat("3.For comparison: compare_gblup_pblup_heritability('trait_name')\n")
}

# Usage examples
cat("=== PBLUP HERITABILITY USAGE EXAMPLES ===\n")
cat("1.With pre-computed A matrix:\n")
cat("   results <- calculate_pblup_heritability_with_Amatrix('your_trait_name')\n\n")

cat("2.From pedigree file:\n")
cat("   results <- calculate_pblup_heritability_from_pedigree('data/pedigree.txt', 'your_trait_name')\n\n")

cat("3.Compare GBLUP vs PBLUP:\n")
cat("   comparison <- compare_gblup_pblup_heritability('your_trait_name')\n\n")

cat("4.Analyze all traits:\n")
cat("   summary <- analyze_all_traits_both_methods()\n\n")

# Run data check
check_pblup_data_availability()

# Method 1: If you have the A matrix already
results_pblup <- calculate_pblup_heritability_with_Amatrix("Is_Dead")

# Method 2: If you only have pedigree file
results_pblup <- calculate_pblup_heritability_from_pedigree("data/pedigree.txt", "Is_Dead")

# Method 3: Compare both methods
comparison <- compare_gblup_pblup_heritability("Is_Dead")

# Method 4: Analyze all traits with both methods
summary_table <- analyze_all_traits_both_methods()
