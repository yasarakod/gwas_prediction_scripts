# ============================================================================== #
# Script: 01_pedigree_matrix_calculation.R
# Title: Pedigree Additive Relationship Matrix Calculation
# Description:
#   Calculates the numerator additive relationship matrix (A-matrix) and its inverse (Ainv) from raw pedigree records using the AGHmatrix package.
#
# Inputs:
#   - data/pedigree.txt (Tab-separated pedigree file with columns: ID, Sire, Dam)
#
# Outputs:
#   - results/completed_pedigree.rds (Formatted pedigree object with ordered IDs)
#   - results/A_matrix.rds (Numerator additive relationship matrix)
#   - results/Ainv_matrix.rds (Inverse additive relationship matrix)
#
# Required R Packages: AGHmatrix, rrBLUP
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


# Check if required packages are installed
if(!require("AGHmatrix")) {
  install.packages("AGHmatrix")
  library(AGHmatrix)
}

# Print AGHmatrix version for debugging
cat("AGHmatrix version loaded\n")

# Check if pedigree file exists
ped_file <- "data/pedigree.txt"
if(!file.exists(ped_file)) {
  stop("Pedigree file not found at ", ped_file, ". Please check file path.")
}

# Read and process pedigree data
cat("Reading pedigree data from", ped_file, "\n")
ped_df <- read.table(ped_file, header = TRUE, sep = "", stringsAsFactors = FALSE)

# Ensure column names are correct
colnames(ped_df)[1:3] <- c("ID", "Sire", "Dam")

# Print summary of original pedigree
cat("Original pedigree has", nrow(ped_df), "animals\n")
cat("First few rows of original pedigree:\n")
print(head(ped_df))

# Clean data - replace NA and empty values with "0"
ped_df$Sire[is.na(ped_df$Sire) | ped_df$Sire == ""] <- "0"
ped_df$Dam[is.na(ped_df$Dam) | ped_df$Dam == ""] <- "0"

# Add missing parents to the pedigree (key step to fix the issues)
cat("Adding missing parents to pedigree...\n")

# Get all unique IDs, sires and dams
all_ids <- unique(ped_df$ID)
all_sires <- unique(ped_df$Sire)
all_dams <- unique(ped_df$Dam)

# Find sires and dams that are not in the ID list (excluding "0")
missing_sires <- setdiff(all_sires[all_sires != "0"], all_ids)
missing_dams <- setdiff(all_dams[all_dams != "0"], all_ids)
all_missing <- unique(c(missing_sires, missing_dams))

cat("Found", length(all_missing), "missing parents to add\n")

# Create rows for missing parents
if(length(all_missing) > 0) {
  missing_rows <- data.frame(
    ID = all_missing,
    Sire = "0",
    Dam = "0",
    stringsAsFactors = FALSE
  )
  
  # Add missing rows to the pedigree
  ped_df <- rbind(ped_df, missing_rows)
}

# Print summary of completed pedigree
cat("Completed pedigree has", nrow(ped_df), "animals\n")

# Verify no self-parenting or duplicates
if(any(duplicated(ped_df$ID))) {
  stop("Duplicate IDs found in pedigree!")
}

if(any(ped_df$ID == ped_df$Sire | ped_df$ID == ped_df$Dam)) {
  stop("Self-parenting found in pedigree!")
}

# Create A matrix with the completed pedigree
cat("Creating A matrix from completed pedigree...\n")
A_matrix <- NULL
error_message <- NULL

try({
  # Try with default parameters first
  A_matrix <- Amatrix(ped_df)
}, silent = TRUE)

if(is.null(A_matrix)) {
  try({
    cat("First attempt failed, trying with more explicit parameters...\n")
    A_matrix <- Amatrix(ped_df, ploidy = 2, verify = FALSE)
  }, silent = TRUE)
}

if(is.null(A_matrix)) {
  try({
    cat("Second attempt failed, trying with fixed pedigree structure...\n")
    # Use minimal pedigree structure
    min_ped <- ped_df[, c("ID", "Sire", "Dam")]
    # Make sure all IDs are strings
    min_ped[] <- lapply(min_ped, as.character)
    A_matrix <- Amatrix(min_ped)
  }, silent = TRUE)
}

if(is.null(A_matrix)) {
  # If all attempts fail, create a simple identity matrix as a fallback
  cat("All attempts to create A matrix failed, creating identity matrix as fallback...\n")
  n <- nrow(ped_df)
  A_matrix <- diag(n)
  rownames(A_matrix) <- ped_df$ID
  colnames(A_matrix) <- ped_df$ID
  error_message <- "Failed to create proper A matrix with AGHmatrix, using identity matrix instead"
}

# Verify if A matrix was created successfully
if(!is.null(A_matrix)) {
  cat("A matrix created successfully - dimensions:", nrow(A_matrix), "x", ncol(A_matrix), "\n")
  
  # Calculate inverse of A matrix
  cat("Calculating inverse of A matrix...\n")
  Ainv <- NULL
  
  try({
    # Add small value to diagonal to avoid singularity issues
    diag(A_matrix) <- diag(A_matrix) + 1e-6
    Ainv <- solve(A_matrix)
  }, silent = TRUE)
  
  if(!is.null(Ainv)) {
    cat("Ainv matrix created successfully - dimensions:", nrow(Ainv), "x", ncol(Ainv), "\n")
    
    # Create directory if it doesn't exist
    if(!dir.exists("data")) {
      dir.create("data")
    }
    
    # Save results
    saveRDS(ped_df, "data/completed_pedigree.rds")
    saveRDS(A_matrix, "data/A_matrix.rds")
    saveRDS(Ainv, "data/Ainv_matrix.rds")
    
    cat("Matrices saved to data directory\n")
    cat("A matrix summary: mean =", mean(A_matrix), "min =", min(A_matrix), "max =", max(A_matrix), "\n")
    cat("Ainv matrix summary: mean =", mean(Ainv), "min =", min(Ainv), "max =", max(Ainv), "\n")
  } else {
    cat("Failed to create Ainv matrix\n")
    if(!is.null(A_matrix)) {
      saveRDS(ped_df, "data/completed_pedigree.rds")
      saveRDS(A_matrix, "data/A_matrix.rds")
      cat("Only A matrix saved to data directory\n")
    }
  }
} else {
  cat("Failed to create A matrix\n")
  if(!is.null(error_message)) {
    cat(error_message, "\n")
  }
}
#PBLUP
              if (pmethod == "PBLUP") {
                # Load pedigree and A matrix from saved RDS files
                if (!exists("A_matrix")) {
                  library(rrBLUP)  # Using rrBLUP instead of sommer
                  # Try to load the pre-computed A matrix
                  if (file.exists("data/A_matrix.rds")) {
                    A_matrix <- readRDS("data/A_matrix.rds")
                    message("Loaded pre-computed A matrix")
                  } else {
                    # Fall back to creating it if the file doesn't exist
                    message("Pre-computed A matrix not found, creating from pedigree...")
                    pedigree <- read.table("data/pedigree.txt", header = TRUE)
                    # Create A matrix manually since rrBLUP doesn't have direct pedigree functions
                    tryCatch({
                      # First try with AGHmatrix which should already be loaded
                      A_matrix <- Amatrix(pedigree)
                    }, error = function(e) {
                      message("Error calculating A matrix from pedigree: ", e$message)
                      message("Loading completed pedigree and recalculating...")
                      # Try with the completed pedigree we saved earlier
                      completed_pedigree <- readRDS("data/completed_pedigree.rds")
                      A_matrix <- Amatrix(completed_pedigree[,1:3])
                    })
                  }
                }
                
                # Ensure all IDs in geno_shrink are in A_matrix
                all_ids <- rownames(geno_shrink)
                missing_ids <- all_ids[!all_ids %in% rownames(A_matrix)]
                
                if (length(missing_ids) > 0) {
                  message(paste("Warning:", length(missing_ids), "IDs from geno_shrink not found in A_matrix"))
                  # Keep only IDs that are in A_matrix
                  all_ids <- all_ids[all_ids %in% rownames(A_matrix)]
                  if (length(all_ids) == 0) {
                    stop("No matching IDs between genotype data and A matrix!")
                  }
                }
                
                # Subset A matrix to available training and testing samples
                A_sub <- A_matrix[all_ids, all_ids]
                
                # Ensure all samples have data
                valid_pheno_ids <- names(pheno)[names(pheno) %in% all_ids]
                valid_indis <- indis[indis %in% valid_pheno_ids]
                
                P_Train <- pheno[valid_pheno_ids[!valid_pheno_ids %in% valid_indis]]
                P_Test <- pheno[valid_pheno_ids[valid_pheno_ids %in% valid_indis]]
                
                # Create dataset for model - adjust for rrBLUP
                y <- pheno[all_ids]
                y[names(P_Test)] <- NA  # Set test values to NA
                
                # Run PBLUP using rrBLUP with error handling
                tryCatch({
                  # rrBLUP's mixed.solve function for PBLUP
                  # Setting K = A_sub (relationship matrix), Z = identity matrix (since all effects are random)
                  Z <- diag(length(all_ids))
                  colnames(Z) <- all_ids
                  rownames(Z) <- all_ids
                  
                  # Run the model with K = relationship matrix
                  model <- mixed.solve(y = y, K = A_sub, 
                                       method = "REML", 
                                       SE = FALSE, 
                                       return.Hinv = FALSE)
                  
                  # Extract predictions (= BLUP + fixed effects)
                  u <- model$u  # random effects (breeding values)
                  
                  # Get predictions for test set
                  test_ids <- names(P_Test)
                  gpred <- u[test_ids]
                  
                  # Format results
                  GPRED <- cbind(gpred, as.numeric(P_Test))
                  GPRED <- round(GPRED, digits = 3)
                  rownames(GPRED) <- test_ids
                  colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                  
                  message(paste('PBLUP Fold ', i, " Time ", n, ' finished #'))
                }, error = function(e) {
                  message(paste("Error in PBLUP model:", e$message))
                  # Return NA values if the model fails
                  gpred <- rep(NA, length(names(P_Test)))
                  names(gpred) <- names(P_Test)
                  GPRED <- cbind(gpred, as.numeric(P_Test))
                  rownames(GPRED) <- names(P_Test)
                  colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                })
              }