# ============================================================================== #
# Script: 05_genomic_prediction_cli_module.R
# Title: Command-Line Modular Genomic Prediction Pipeline
# Description:
#   Provides a flexible command-line interface using argparse to run cross-validation prediction workflows with configurable traits, folds, iterations, and model parameters.
#
# Inputs:
#   - data/pheno (Phenotype file)
#   - data/M (Genotype file)
#   - Command-line flags (--dir, --pheno_file, --geno_file, --fold, --iterations, --trait_name, --marker_size)
#
# Outputs:
#   - results/<Method>_Folds_<fold>_Times_<iter>_<trait>_<marker_size>.txt (Prediction results)
#
# Required R Packages: argparse, rrBLUP, BGLR, glmnet, randomForest, e1071, brnn, AGHmatrix
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




# ---- Argument Handling (Combined Approach) ----
library(argparse)

parser <- ArgumentParser(description = "CV Analysis Module")

# Define Arguments (Keep this section unchanged)
parser$add_argument("--dir", default = "data", required = TRUE,help = "Data directory")
parser$add_argument("--pheno_file", default = "pheno", help = "Phenotype file name")
parser$add_argument("--geno_file", default = "M", help = "Genotype file name")
parser$add_argument("--fold", type = "integer", default = 3, help = "Number of folds")
parser$add_argument("--iter", type = "integer", default = 3, help = "Number of iterations")
parser$add_argument("--trait_names", nargs = "+",required = TRUE, help = "Trait names (space-separated)")
parser$add_argument("--pmethods", nargs = "+", default = c("PBLUP","GBLUP","BA","BB","BC","BL","BRR","EN","RR","RF"), help = "Prediction methods (space-separated)")
parser$add_argument("--Fixed_effect", nargs = "+", default = c("NULL"), help = "Fixed effects (space-separated)")
parser$add_argument("--pop_sizes", nargs = "+", default = c("NULL"), help = "Population sizes (space-separated)")
parser$add_argument("--marker_sizes", nargs = "+", default = c("NULL"), help = "Marker sizes (space-separated)")

# Parse Arguments
args <- parser$parse_args()

data_path = paste0(args$dir, "/")

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


import_packages<-function(){
  #define packages to install
  packages <- c("data.table","snpReady","SNPRelate","rrBLUP","BGLR","dplyr","glmnet","randomForest","e1071","brnn","pedigree"
                ,"reshape2","regress","MASS")
  
  #if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  
  #BiocManager::install("SNPRelate")
  #BiocManager::install("snpReady")
  #install all packages that are not already installed
  install.packages(setdiff(packages, rownames(installed.packages())))
  
  invisible(lapply(packages, library, character.only=TRUE))
}

import_packages()

run_CV<- function(pheno_file, geno_file, fold,iter,trait_names, pmethods, Fixed_effect,pop_sizes,marker_sizes,data_path){
  
  # --- Load and Preprocess Data ---
  pheno_data <- readRDS(paste0(data_path,pheno_file))
  geno_data <- readRDS(paste0(data_path,geno_file))
  
  # Initialize a list to store accuracy results
  accuracy_results <- list()
  # --- Main Loop for Trait, Method, Sizes ---
  for (popsize in pop_sizes){
    for (markersize in marker_sizes) {
      for (pmethod in pmethods) {
        for (trait in trait_names){
          # Find the correct column for the current trait
          pheno.col <- which(colnames(pheno_data) == trait)
          message(paste0("Running genomic prediction for trait '", trait, 
                         "' with method '", pmethod, "' and marker size: ", markersize," and pop size: ",popsize))
          ped_formated<-geno_data
          phe<-pheno_data
          pheno <- setNames(as.numeric(phe[, pheno.col]), phe[, 1])
          
          # Find common samples with complete data
          common_samples=intersect(rownames(ped_formated),names(pheno))
          
          if (length(common_samples) == 0) {
            stop(paste0("No common samples with complete data for trait '", trait, "'."))
          }
          
          # Sampling population (optional) and subsetting data
          if (popsize != "NULL") {
            ped_formated <- ped_formated[sample(common_samples, popsize), , drop = FALSE]
            pheno <- pheno[rownames(ped_formated)]
          }
          
          # Informative messages
          message("Number of common samples between geno and pheno: ", nrow(ped_formated))
          message("Number of markers after filtering: ", ncol(ped_formated))
          
          # Function to randomly sample markers
          geno_shrink <- ped_formated[, if (!is.null(markersize) && markersize < ncol(ped_formated)) 
            sample(ncol(ped_formated), markersize) 
            else 1:ncol(ped_formated), drop = FALSE]
          
          message("New genotypic data dimension: ", dim(geno_shrink)[1], " x ", dim(geno_shrink)[2])
          
          # Handle fixed effects (if any)
          FIXED <- if (Fixed_effect[1] != "NULL") {
            df <- phe[, Fixed_effect, drop = FALSE]
            rownames(df) <- phe$Sample_ID  # Use rownames for consistency
            df[rownames(ped_formated), , drop = FALSE] 
          } else {
            NULL
          }
          
          # Save filtered data
          saveRDS(pheno, "pheno.rds")
          saveRDS(geno_shrink, "geno_shrink.rds")
          
          df = list()
          for (n in 1:iter){ 
            df_total = data.frame()
            folds <- list() # flexible object for storing folds
            fold.size <- length(pheno)/fold
            remain <- names(pheno) # all obs are in
            
            for (i in 1:fold){
              select <- sample(remain, fold.size, replace = FALSE)
              folds[[i]] <- select # store indices
              
              if (i == fold){
                folds[[i]] <- remain
              }
              
              remain <- setdiff(remain, select)
              remain
              
              indis <- folds[[i]] #unpack into a vector
              P_Train<-data.frame(pheno[!(names(pheno) %in% indis) ] )
              P_Test <- data.frame(pheno[names(pheno) %in% indis ] )
              
              G_Train<-as.matrix(geno_shrink[rownames(geno_shrink) %in% rownames(P_Train), ] )
              G_pred<-as.matrix(geno_shrink[!rownames(geno_shrink) %in% rownames(P_Train), ] )
              
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
              
              
              #GBLUP
              if (pmethod == "GBLUP"){
                
                if(Fixed_effect[1]!="NULL"){
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  test<-mixed.solve(y=P_Train[,1],Z=G_Train,X = as.matrix(cbind(matrix(1, nrow(G_Train)),as.numeric(FixedTrain))))
                  MODEL<- as.vector(G_pred %*% as.matrix(test$u))+ as.matrix(matrix(1, nrow(G_pred))) %*% test$beta
                }
                else{
                  test<-mixed.solve(y=P_Train[,1],Z=G_Train,X = matrix(1, nrow(G_Train)))
                  MODEL<- as.vector(G_pred %*% as.matrix(test$u))+ as.matrix(matrix(1, nrow(G_pred))) %*% test$beta
                }
                
                gpred = data.frame(MODEL)[,1]
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("GBLUP_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('GBLUP Fold ',i," Time ",n,' finished #'))
              }
              
              #BayesA(Scaled-t prior)
              if (pmethod == "BA"){
                
                GENO = rbind(G_Train,G_pred)
                colnames(P_Train)<-"V"
                colnames(P_Test)<-"V"
                y<- rbind(P_Train,P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa<-setNames(y$V,rownames(y))
                nIter=1500;
                burnIn=500;
                #thin=5;
                #saveAt ='';
                #S0=NULL;
                #weights=NULL;
                #R2=0.5;
                
                if(Fixed_effect[1]!="NULL"){
                  
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred<-as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain,FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesA'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                else{
                  ETA<-list(list(X=GENO,model='BayesA'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                
                gpred = MODEL$yHat
                gpredSD=MODEL$SD.yHat
                VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
                SDU=MODEL$SD.yHat
                
                CD=sqrt(1-(SDU^2/VARU))
                GPRED=cbind(data.frame(data.frame(gpred)[match(rownames(P_Test),rownames(data.frame(gpred))),]),data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BA_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BA Fold ',i," Time ",n,' finished #'))
                
              }
              
              if (pmethod == "BB"){
                
                GENO = rbind(G_Train,G_pred)
                colnames(P_Train)<-"V"
                colnames(P_Test)<-"V"
                y<- rbind(P_Train,P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa<-setNames(y$V,rownames(y))
                nIter=1500;
                burnIn=500;
                #thin=5;
                #saveAt ='';
                #S0=NULL;
                #weights=NULL;
                #R2=0.5;
                
                if(Fixed_effect[1]!="NULL"){
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred<-as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain,FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesB'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                else{
                  options(warn=-1)
                  ETA<-list(list(X=GENO,model='BayesB'))
                  MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                
                gpred = MODEL$yHat
                gpredSD=MODEL$SD.yHat
                VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
                SDU=MODEL$SD.yHat
                
                CD=sqrt(1-(SDU^2/VARU))
                GPRED=cbind(data.frame(data.frame(gpred)[match(rownames(P_Test),rownames(data.frame(gpred))),]),data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BB_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BB Fold ',i," Time ",n,' finished #'))
                
              }
              
              if (pmethod == "BC"){
                
                GENO = rbind(G_Train,G_pred)
                colnames(P_Train)<-"V"
                colnames(P_Test)<-"V"
                y<- rbind(P_Train,P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa<-setNames(y$V,rownames(y))
                nIter=1500;
                burnIn=500;
                #thin=5;
                #saveAt ='';
                #S0=NULL;
                #weights=NULL;
                #R2=0.5;
                
                if(Fixed_effect[1]!="NULL"){
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred<-as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain,FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BayesC'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                else{
                  options(warn=-1)
                  ETA<-list(list(X=GENO,model='BayesC'))
                  MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                
                gpred = MODEL$yHat
                gpredSD=MODEL$SD.yHat
                VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
                SDU=MODEL$SD.yHat
                
                CD=sqrt(1-(SDU^2/VARU))
                GPRED=cbind(data.frame(data.frame(gpred)[match(rownames(P_Test),rownames(data.frame(gpred))),]),data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BC_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BC Fold ',i," Time ",n,' finished #'))
                
              }
              
              if (pmethod == "BL"){
                
                GENO = rbind(G_Train,G_pred)
                colnames(P_Train)<-"V"
                colnames(P_Test)<-"V"
                y<- rbind(P_Train,P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa<-setNames(y$V,rownames(y))
                nIter=1500;
                burnIn=500;
                #thin=5;
                #saveAt ='';
                #S0=NULL;
                #weights=NULL;
                #R2=0.5;
                
                if(Fixed_effect[1]!="NULL"){
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred<-as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain,FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BL'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                else{
                  options(warn=-1)
                  ETA<-list(list(X=GENO,model='BL'))
                  MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                
                gpred = MODEL$yHat
                gpredSD=MODEL$SD.yHat
                VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
                SDU=MODEL$SD.yHat
                
                CD=sqrt(1-(SDU^2/VARU))
                GPRED=cbind(data.frame(data.frame(gpred)[match(rownames(P_Test),rownames(data.frame(gpred))),]),data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BL_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BL Fold ',i," Time ",n,' finished #'))
                
              }
              
              #BayesRR
              if (pmethod == "BRR"){
                
                GENO = rbind(G_Train,G_pred)
                colnames(P_Train)<-"V"
                colnames(P_Test)<-"V"
                y<- rbind(P_Train,P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa<-setNames(y$V,rownames(y))
                nIter=1500;
                burnIn=500;
                
                
                if(Fixed_effect[1]!="NULL"){
                  FixedTrain<-as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred<-as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain,FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2<-list(list(X=FIX,model="FIXED"),list(X=GENO,model='BRR'))
                  options(warn=-1)
                  MODEL=BGLR(y=yNa,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                else{
                  options(warn=-1)
                  ETA<-list(list(X=GENO,model='BRR'))
                  MODEL=BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose = F)
                }
                
                gpred = MODEL$yHat
                gpredSD=MODEL$SD.yHat
                VARU=var(gpred)   # See: BGLR Genomics.cimmyt.org/BGLR-extdoc.pdf
                SDU=MODEL$SD.yHat
                
                CD=sqrt(1-(SDU^2/VARU))
                GPRED=cbind(data.frame(data.frame(gpred)[match(rownames(P_Test),rownames(data.frame(gpred))),]),data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BRR_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BRR Fold ',i," Time ",n,' finished #'))
              }
              
              if (pmethod == "EN"){
                
                cv.fit <- cv.glmnet(G_Train,P_Train[,1],family="gaussian",alpha=0.5,nfolds=10) #ElasticNet penalty with top results
                lambda_min <- cv.fit$lambda.min # making the best prediction
                gpred <- predict(cv.fit,newx=G_pred,s=c(lambda_min))
                #GPRED <- round(gpred,digits=3)
                gpredSD=gpred
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("EN_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('EN Fold ',i," Time ",n,' finished #'))
              }
              
              if (pmethod == "BRNN"){
                
                cv.fit <- brnn(G_Train,P_Train[,1],neurons=2,epochs=10, verbose=T) # neural networks with 2 neurons and 50 iterations
                gpred <- predict(cv.fit,newdata=G_pred)
                
                gpredSD=gpred
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("BRNN_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('BRNN Fold ',i," Time ",n,' finished #'))
              }
              
              if (pmethod == "RR"){
                
                cv.fit <- cv.glmnet(G_Train,P_Train[,1],alpha=0) #Ridge penalty with top results
                
                lambda_min <- cv.fit$lambda.min # making the best prediction
                
                gpred <- predict(cv.fit,newx=G_pred,s=c(lambda_min))
                
                gpredSD=gpred
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("RR_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('RR Fold ',i," Time ",n,' finished #'))
              }
              
              if (pmethod == "RF"){
                
                gpred = randomForest(G_Train, P_Train[,1], xtest=G_pred)$test$predicted;
                
                gpredSD=gpred
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("RF_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('RF Fold ',i," Time ",n,' finished #'))
              }
              
              if (pmethod == "SVM"){
                
                model <- svm(G_Train,P_Train,method="nu-regression",kernel="radial",cost=10,gamma=0.001)
                gpred <- predict(model,G_pred)
                
                gpredSD=gpred
                gpredSD=NA
                CD <- gpred
                CD <- NA
                GPRED=cbind(gpred,data.frame(P_Test))
                GPRED=round(GPRED,digits=3)
                rownames(GPRED)<-rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted',n),paste0('Actual',n))
                #write.table(GPRED,paste("SVM_Fold_",i,"_Time_",n,"_",sname,"_",paste(Fixed_effect,collapse = "_"),".txt",sep = ""),quote = F, sep="\t")
                message(paste('SVM Fold ',i," Time ",n,' finished #'))
              }
              
              GPRED<-as.matrix(GPRED)
              df_total<-rbind(df_total,GPRED)
              
            }
            
            df[[n]]<-df_total
            
          }
          
          df_t<-data.frame()
          for (r in 1:iter){
            if(r==1){
              df_t<-data.frame(df[[r]])
            }
            else{
              df_m<-data.frame(df[[r]])
              df_mt<-df_m[match(rownames(df_t),rownames(df_m)),]
              df_t<-cbind(df_t,df_mt)
            }
          }
          predicted_cols <- grep("Predicted", names(df_t))  # Find all columns starting with "Predicted"
          actual_cols <- grep("Actual", names(df_t))       # Find all columns starting with "Actual"
          
          # Calculate and add mean columns
          df_t$PredictedMean <- rowMeans(df_t[, predicted_cols])
          df_t$ActualMean <- rowMeans(df_t[, actual_cols])
          
          # Round the mean columns
          df_t[, c("PredictedMean", "ActualMean")] <- round(df_t[, c("PredictedMean", "ActualMean")], digits = 3)
          df1<-df_t
          write.table(df1,paste(data_path,pmethod,"_Folds_",i,"_Times_",n,"_",trait,"_",paste(Fixed_effect,collapse = "_"),"_",popsize,"_",markersize,".txt",sep = ""),quote = F, sep="\t")
          print(cor(df1$PredictedMean,df1$ActualMean))
          accuracy_results[[paste(trait, pmethod, popsize, markersize)]] <- cor(df1$PredictedMean,df1$ActualMean)
          
        } #end of trait names
      } #end of pmethods
    } #end of marker sizes
  } #end of pop sizes
  
  Accuracy = unlist(accuracy_results)
  write.table(Accuracy,paste0(data_path,"accuracy_results_M",Sys.Date(),".txt"),sep = "\t",quote = F, col.names = F)
}

run_CV(pheno_file=args$pheno_file,
       geno_file=args$geno_file,
       fold = args$fold,
       iter = args$iter,
       trait_names = args$trait_names,
       pmethods= args$pmethods,
       Fixed_effect = args$Fixed_effect,
       pop_sizes = args$pop_sizes,
       marker_sizes = args$marker_sizes,
       data_path=data_path
       )
