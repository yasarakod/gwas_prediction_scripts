# ============================================================================== #
# Script: 06_weighted_gblup_prediction.R
# Title: Iterative Weighted GBLUP (wGBLUP) Pipeline
# Description:
#   Implements iterative Weighted GBLUP where SNP marker weights are calculated from estimated marker effects to construct trait-specific weighted genomic relationship matrices (G-matrix).
#
# Inputs:
#   - data/M (Genotype dosage matrix)
#   - data/pheno (Phenotype dataset)
#   - data/wald_*.txt (GWAS Wald test summary statistics for initial weights)
#
# Outputs:
#   - results/WGBLUP_Folds_*_Times_*.txt (Cross-validation prediction accuracy logs)
#   - results/marker_weights_<trait>_<marker_size>.rds (Calculated SNP variance weights)
#
# Required R Packages: data.table, snpReady, SNPRelate, rrBLUP, BGLR, dplyr, glmnet, randomForest, e1071, brnn, pedigree, reshape2, regress, MASS
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

################################################################################
#process ped data
################################################################################
process_ped_data <- function(input_file, output_file) {
  # Read the data from the PED file
  ped <- read.table(file = input_file, sep = "\t")
  
  # Define a vector of regex patterns and replacements for cleaner code
  regex_patterns <- c(
    ".*dF0(.*)\\.CEL"
  )
  
  replacements <- c(
    "dF0\\1"
  )
  
  # Loop through patterns and apply replacements
  for (i in seq_along(regex_patterns)) {
    ped$V1 <- sub(regex_patterns[i], replacements[i], ped$V1)
  }
  
  # Write the processed data back to the output files
  write.table(ped$V1, "raw_file_names_original.txt", quote = FALSE)
  write.table(ped, output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# Usage
input_file <- "DGS_dF0_H.ped"
output_file <- "DGS_dF0_H.ped"
process_ped_data(input_file, output_file)

################################################################################
#process_map_data
################################################################################
process_map_data <- function(input_file, output_file) {
  # Read the data from the map file
  map <- read.table(input_file, sep = "\t")
  
  chr <- c("CM007757.1", "CM007758.1", "CM007759.1", "CM007760.1", "CM007761.1", "CM007762.1",
           "CM007763.1", "CM007764.1", "CM007765.1", "CM007766.1", "CM007767.1", "CM007768.1",
           "CM007769.1", "CM007770.1", "CM007771.1", "CM007772.1", "CM007773.1", "CM007774.1",
           "CM007775.1", "CM007776.1", "CM007777.1", "CM007778.1", "CM007779.1", "CM007780.1")
  
  cont <- unique(map$V1)
  
  # Replace values from 'chr' vector with corresponding numeric values
  for (i in 1:length(chr)) {
    map$V1 <- gsub(pattern = chr[i], replacement = i, map$V1)
  }
  
  # Replace 'cont' values with "contigX" where X is the index
  for (i in 1:length(cont)) {
    if (nchar(cont[i]) < 3) {}
    else {map$V1 <- gsub(cont[i], paste0("contig", i), map$V1)
    }
  }
  
  # Write the processed data back to the output file
  write.table(map$V1, "contigs.txt", sep = "\t", quote = FALSE)
  write.table(map, output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# usage
input_file <- "DGS_dF0_H.map"
output_file <- "DGS_dF0_H.map"
process_map_data(input_file, output_file)

################################################################################
#QC filtering
################################################################################

process_plink_commands <- function() {
  # Convert to binary format and perform initial QC
  shell("plink.exe --aec --autosome-num 24 --file DGS_dF0_H --make-bed --out DGS_dF0_binary --no-fid --no-parents --no-sex --allow-no-sex --no-pheno")
  
  # # Read and filter the . fam file
  # fam <- read.table("DGS_dF0_binary.fam", sep = " ")
  # fam1 <- fam[grep("^dF0", fam$V1, perl = TRUE), ]
  # write.table(fam1, "samples.fam", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Perform further filtering with PLINK
  shell("plink.exe --aec --autosome-num 24 --bfile DGS_dF0_binary --keep selected.txt --make-bed --out DGS_dF0_filtered")
  
  # Additional quality control
  shell("plink.exe --aec --autosome-num 24 --bfile DGS_dF0_filtered --make-bed --geno 0.1 --maf 0.05 --out DGS_dF0_qc")
}

# Usage
process_plink_commands()

# setwd() configured via project environment
library(data.table)
library(snpReady)
library(SNPRelate)
library(base)

args <- c("data","DGS_dF0_qc","selected.txt","phenotypes.csv")

data_path = paste0(args[1], "/")
base_filename <- args[2]
selection<- args[3]
phenotype<-args[4]


preprocess_genotype <- function(infile, selection,gwas_snp) {
  #detach("package: dplyr", unload = TRUE)
  base:: system(paste0("plink.exe --allow-extra-chr --keep ",selection," --bfile ",infile," --make-bed --out ",data_path,"tempgeno"))
  base::system(paste0("plink.exe --allow-extra-chr --bfile ",data_path,"tempgeno --recode A --out ",data_path,"tempgeno"))
  base::system(paste0("plink.exe --allow-extra-chr --bfile ",data_path,"tempgeno --recode -tab --out ",data_path,"tempgeno"))
  
  ped<-dataframe(fread(paste0(data_path,"tempgeno.raw"),header = T)[,-c(1,3: 6)])
  rownames(ped)<-ped$IID
  ped<-ped[,-c(1)]
  colnames(ped)<-gsub("_[ATCG]$","",colnames(ped))
  ped<-as.matrix(ped)
  
  hapmap <- read.table(paste0(data_path,"tempgeno.map"))[,c(2, 1, 4)]
  colnames(hapmap)<-c("rs","chrom","pos")
  hapmap$rs<-gsub("-",".",hapmap$rs)
  hapmap$chrom <- as.factor(hapmap$chrom)
  hapmap$pos <- as.numeric(hapmap$pos)
  return(list(ped = ped, hapmap = hapmap))  # Return results in a list
}

perform_marker_qc <- function(ped, hapmap) {
  QC <- raw.data(data = ped, frame = "wide",hapmap = hapmap,
                 sweep.sample = 1, maf = 0.05, call.rate = 0.90, 
                 base = FALSE, imput = TRUE, imput.type = "wright", 
                 outfile = "012", plot = F)
  
  # Ensure markers and map are in the same order (this check is important)
  stopifnot(identical(as.character(QC$Hapmap$rs), colnames(QC$M.clean)))
  
  return(list(M = QC$M.clean, hapmap = QC$Hapmap))  # Return results in a list
}

prune_markers_by_ld <- function(M, hapmap) {
  
  # creating the GDS file
  snpgdsCreateGeno(gds.fn = "toy.gds", genmat = M,sample.id = rownames(M), 
                   snp.id = as.integer(1:dim(hapmap)[1]),snp.rs.id = hapmap$rs,                 
                   snp.chromosome = as.integer(hapmap$chrom),snp.position = as.integer(hapmap$pos),     
                   snpfirstdim = FALSE)                  
  
  # Loading gds
  genofile <- snpgdsOpen(filename = "toy.gds")
  
  # prune markers by MAF 0.10, CR 0.90 e LD 0.99 (r2)
  snps_pruned <- snpgdsLDpruning(gdsobj = genofile,remove.monosnp = TRUE,     
                                 maf = 0.05,missing.rate = 0.90,method = "corr",
                                 slide.max.bp = 100000,ld.threshold = 0.98,start.pos = "first",
                                 autosome.only=F)       
  
  # get SNP ids
  snps_pruned <- unlist(unname(snps_pruned))
  
  # removing SNPs prunned by LD from M matrix
  M <- M[, as.numeric(snps_pruned)]
  
  # correcting the map object
  hapmap <- hapmap[as.numeric(snps_pruned), ] # removing SNPs pruned by LD
  
  # verifying if all the markers into M matrix are in the same order into map
  identical(as.character(hapmap$rs), colnames(M))
  
  # close GDS
  snpgdsClose(genofile)
  
  # saving files
  saveRDS(hapmap, paste0(data_path,"hapmap"))
  saveRDS(M, paste0(data_path,"M"))
  return(list(M = M, hapmap = hapmap))  # Return results in a list
}

match_phenotype <- function(M, pheno_file) {
  # Check if phenotype file exists
  if (!file.exists(paste0(data_path,pheno_file))) {
    stop("Phenotype file not found:  ", pheno_file)
  }
  
  pheno <- read.csv(paste0(data_path,pheno_file), header = TRUE)
  
  # Check for matching sample IDs
  if (! all(rownames(M) %in% pheno$Sample_ID)) {
    stop("Not all sample IDs in genotype data found in phenotype data.")
  }
  
  pheno_matched <- pheno[match(rownames(M), pheno$Sample_ID),]
  saveRDS(pheno_matched, paste0(data_path,"pheno"))
  return(pheno_matched)
}

geno_data <- preprocess_genotype(paste0(data_path,base_filename), paste0(data_path,selection))
geno_data <- perform_marker_qc(geno_data$ped, geno_data$hapmap)
geno_data <- prune_markers_by_ld(geno_data$M, geno_data$hapmap)
pheno_matched <- match_phenotype(geno_data$M, phenotype)

################################################################################
################################################################################
import_packages<-function(){
  #define packages to install
  packages <- c("data.table","snpReady","SNPRelate","rrBLUP","BGLR","dplyr","glmnet","randomForest","e1071","brnn","pedigree"
                ,"reshape2","regress","MASS")
  
  install.packages(setdiff(packages, rownames(installed.packages())))
  
  invisible(lapply(packages, library, character.only=TRUE))
}

import_packages()

################################################################################
# NEW FUNCTION:  Compute marker weights from GWAS p-values
################################################################################
compute_marker_weights <- function(snp_pvalues, method = "inverse_pvalue", alpha = 1) {
  # Method options: 
  # - "inverse_pvalue": weights = 1/p or -log10(p)
  # - "bonferroni": weights based on significance threshold
  # - "sqrt_inverse":  weights = sqrt(1/p)
  # - "power":  weights = (1/p)^alpha
  
  if (method == "inverse_pvalue") {
    # Use -log10(p) as weights (handles very small p-values better)
    weights <- -log10(snp_pvalues$p + 1e-300) # add small value to avoid log(0)
  } else if (method == "sqrt_inverse") {
    weights <- sqrt(1 / (snp_pvalues$p + 1e-300))
  } else if (method == "power") {
    weights <- (1 / (snp_pvalues$p + 1e-300))^alpha
  } else if (method == "bonferroni") {
    # Binary weights based on Bonferroni threshold
    threshold <- 0.05 / nrow(snp_pvalues)
    weights <- ifelse(snp_pvalues$p < threshold, 1, 0.1) # significant SNPs get weight 1, others 0.1
  } else {
    stop("Unknown weighting method")
  }
  
  # Normalize weights to have mean of 1 (optional but recommended)
  weights <- weights / mean(weights)
  
  return(data.frame(id = snp_pvalues$id, weight = weights))
}

################################################################################
# MODIFIED FUNCTION: Weighted GBLUP using weighted relationship matrix
################################################################################
weighted_gblup <- function(G_Train, G_pred, P_Train, marker_weights, FIXED = NULL) {
  # G_Train: genotype matrix for training (n x m)
  # G_pred: genotype matrix for prediction (n_pred x m)
  # P_Train: phenotype vector for training
  # marker_weights: vector of weights for each marker (length m)
  # FIXED: optional fixed effects matrix
  
  # Ensure weights match markers
  if (length(marker_weights) != ncol(G_Train)) {
    stop("Number of weights must match number of markers")
  }
  
  # Apply weights to genotype matrices
  G_Train_weighted <- sweep(G_Train, 2, sqrt(marker_weights), "*")
  G_pred_weighted <- sweep(G_pred, 2, sqrt(marker_weights), "*")
  
  # Compute weighted genomic relationship matrix (GRM)
  # G = ZDZ'/sum(2*p*(1-p)*d) where D is diagonal matrix of weights
  p_freq <- colMeans(G_Train) / 2
  norm_factor <- sum(2 * p_freq * (1 - p_freq) * marker_weights)
  
  # Use weighted genotypes for mixed.solve
  if (! is.null(FIXED)) {
    fit <- mixed.solve(y = P_Train[,1], 
                       Z = G_Train_weighted, 
                       X = as.matrix(cbind(matrix(1, nrow(G_Train_weighted)), as.numeric(FIXED))))
    predictions <- as.vector(G_pred_weighted %*% as.matrix(fit$u)) + 
      as.matrix(matrix(1, nrow(G_pred_weighted))) %*% fit$beta
  } else {
    fit <- mixed.solve(y = P_Train[,1], 
                       Z = G_Train_weighted, 
                       X = matrix(1, nrow(G_Train_weighted)))
    predictions <- as.vector(G_pred_weighted %*% as.matrix(fit$u)) + 
      as.matrix(matrix(1, nrow(G_pred_weighted))) %*% fit$beta
  }
  
  return(predictions)
}

run_CV<- function(pheno_file, geno_file, fold, iter, trait_names, pmethods, Fixed_effect, 
                  pop_sizes, marker_sizes, data_path, weight_method = "inverse_pvalue", weight_alpha = 1){
  
  # --- Load and Preprocess Data ---
  pheno_data <- readRDS(paste0(data_path, pheno_file))
  geno_data <- readRDS(paste0(data_path, geno_file))
  
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
                         "' with method '", pmethod, "' and marker size:  ", markersize, 
                         " and pop size: ", popsize))
          ped_formated <- geno_data
          
          phe <- pheno_data
          pheno <- setNames(as.numeric(phe[, pheno.col]), phe[, 1])
          
          # Find common samples with complete data
          common_samples <- intersect(rownames(ped_formated), names(pheno))
          
          if (length(common_samples) == 0) {
            stop(paste0("No common samples with complete data for trait '", trait, "'. "))
          }
          
          # Sampling population (optional) and subsetting data
          if (popsize != "NULL") {
            ped_formated <- ped_formated[sample(common_samples, popsize), , drop = FALSE]
            pheno <- pheno[rownames(ped_formated)]
          }
          
          # Informative messages
          message("Number of common samples between geno and pheno: ", nrow(ped_formated))
          message("Number of markers after filtering: ", ncol(ped_formated))
          
          # Load GWAS p-values for weighting
          snp_data <- read.table(paste0(data_path, "wald_", trait, ".txt"), header = TRUE)
          
          # Sort by p-value (ascending)
          snp_data_sorted <- snp_data[order(snp_data$p), ]
          
          # Select top SNPs
          top_snps <- list(
            top_50 = snp_data_sorted$id[1:50],
            top_60 = snp_data_sorted$id[1:60],
            top_70 = snp_data_sorted$id[1:70],
            top_80 = snp_data_sorted$id[1:80],
            top_90 = snp_data_sorted$id[1:90],
            top_100 = snp_data_sorted$id[1:100],
            top_500 = snp_data_sorted$id[1:500],
            top_1000 = snp_data_sorted$id[1:1000],
            top_5000 = snp_data_sorted$id[1:5000],
            top_10000 = snp_data_sorted$id[1:10000],
            top_20000 = snp_data_sorted$id[1:20000],
            top_30000 = snp_data_sorted$id[1:30000],
            top_40000 = snp_data_sorted$id[1:40000],
            top_50000 = snp_data_sorted$id[1:50000],
            top_ALL = snp_data_sorted$id[1:55158]
          )
          
          # Save each set of top SNPs
          for (name in names(top_snps)) {
            write.table(top_snps[[name]], 
                        file = paste0(data_path, name, "_snps.txt"), 
                        quote = FALSE, 
                        row.names = FALSE, 
                        col.names = FALSE)
          }
          
          # Function to selectively sample markers or use top SNPs from files
          geno_data_shrink <- function(ped_formated, markersize, data_path) {
            if (is.character(markersize) && grepl("^top_", markersize)) {
              # Read top SNPs from file
              top_snps <- read.table(paste0(data_path, markersize, "_snps.txt"), header = FALSE)$V1
              top_snps <- gsub("-", ".", top_snps)
              
              # Subset based on top SNPs
              geno_data_shrink <- ped_formated[, colnames(ped_formated) %in% top_snps, drop = FALSE]
              
              message("Using top SNPs for marker size:  ", markersize)
            } else if (markersize == "ALL") {
              # Use all markers
              geno_data_shrink <- ped_formated
              
              message("Using all markers.")
            } else {
              # Randomly sample markers if markersize is numeric
              geno_data_shrink <- ped_formated[, if (! is.null(markersize) && markersize < ncol(ped_formated)) 
                sample(ncol(ped_formated), markersize) 
                else 1:ncol(ped_formated), drop = FALSE] 
            }
            
            return(geno_data_shrink)
          }
          
          # Call the geno_shrink function and assign the output to a new variable
          geno_shrink <- geno_data_shrink(ped_formated, markersize, data_path)
          
          # Use the new variable to print the dimensions
          message("New genotypic data dimension: ", dim(geno_shrink)[1], " x ", dim(geno_shrink)[2])
          
          # Compute marker weights based on GWAS p-values
          # Match SNP IDs between geno_shrink and snp_data
          snp_data$id <- gsub("-", ".", snp_data$id)
          matched_snps <- snp_data[snp_data$id %in% colnames(geno_shrink), ]
          matched_snps <- matched_snps[match(colnames(geno_shrink), matched_snps$id), ]
          
          # Compute weights
          marker_weight_df <- compute_marker_weights(matched_snps, method = weight_method, alpha = weight_alpha)
          marker_weights <- marker_weight_df$weight
          names(marker_weights) <- marker_weight_df$id
          
          # Ensure weights are in same order as markers in geno_shrink
          marker_weights <- marker_weights[colnames(geno_shrink)]
          
          message("Marker weights computed using method: ", weight_method)
          message("Weight range: ", round(min(marker_weights), 3), " to ", round(max(marker_weights), 3))
          
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
          saveRDS(marker_weights, paste0(data_path, "marker_weights_", trait, "_", markersize, ". rds"))
          
          df <- list()
          for (n in 1:iter){ 
            df_total <- data.frame()
            folds <- list() # flexible object for storing folds
            fold.size <- length(pheno) / fold
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
              P_Train <- data.frame(pheno[!(names(pheno) %in% indis)])
              P_Test <- data.frame(pheno[names(pheno) %in% indis])
              
              G_Train <- as.matrix(geno_shrink[rownames(geno_shrink) %in% rownames(P_Train), ])
              G_pred <- as.matrix(geno_shrink[! rownames(geno_shrink) %in% rownames(P_Train), ])
              
              # Get weights for current markers
              current_weights <- marker_weights[colnames(G_Train)]
              
              #GBLUP
              if (pmethod == "GBLUP"){
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  test <- mixed.solve(y = P_Train[,1], Z = G_Train, X = as.matrix(cbind(matrix(1, nrow(G_Train)), as.numeric(FixedTrain))))
                  MODEL <- as.vector(G_pred %*% as.matrix(test$u)) + as.matrix(matrix(1, nrow(G_pred))) %*% test$beta
                } else {
                  test <- mixed.solve(y = P_Train[,1], Z = G_Train, X = matrix(1, nrow(G_Train)))
                  MODEL <- as.vector(G_pred %*% as.matrix(test$u)) + as.matrix(matrix(1, nrow(G_pred))) %*% test$beta
                }
                
                gpred <- data.frame(MODEL)[,1]
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('GBLUP Fold ', i, " Time ", n, ' finished #'))
              }
              
              # WEIGHTED GBLUP (NEW METHOD)
              if (pmethod == "WGBLUP"){
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  MODEL <- weighted_gblup(G_Train, G_pred, P_Train, current_weights, FIXED = FixedTrain)
                } else {
                  MODEL <- weighted_gblup(G_Train, G_pred, P_Train, current_weights, FIXED = NULL)
                }
                
                gpred <- as.vector(MODEL)
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('WGBLUP Fold ', i, " Time ", n, ' finished #'))
              }
              
              #BayesA(Scaled-t prior)
              if (pmethod == "BA"){
                
                GENO <- rbind(G_Train, G_pred)
                colnames(P_Train) <- "V"
                colnames(P_Test) <- "V"
                y <- rbind(P_Train, P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa <- setNames(y$V, rownames(y))
                nIter <- 1500
                burnIn <- 500
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred <- as.matrix(FIXED[! rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain, FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2 <- list(list(X = FIX, model = "FIXED"), list(X = GENO, model = 'BayesA'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA2, nIter = nIter, burnIn = burnIn, verbose = F)
                } else {
                  ETA <- list(list(X = GENO, model = 'BayesA'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = F)
                }
                
                gpred <- MODEL$yHat
                gpredSD <- MODEL$SD.yHat
                VARU <- var(gpred)
                SDU <- MODEL$SD.yHat
                
                CD <- sqrt(1 - (SDU^2 / VARU))
                GPRED <- cbind(data.frame(data.frame(gpred)[match(rownames(P_Test), rownames(data.frame(gpred))),]), data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BA Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "BB"){
                
                GENO <- rbind(G_Train, G_pred)
                colnames(P_Train) <- "V"
                colnames(P_Test) <- "V"
                y <- rbind(P_Train, P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa <- setNames(y$V, rownames(y))
                nIter <- 1500
                burnIn <- 500
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred <- as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain, FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2 <- list(list(X = FIX, model = "FIXED"), list(X = GENO, model = 'BayesB'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA2, nIter = nIter, burnIn = burnIn, verbose = F)
                } else {
                  options(warn = -1)
                  ETA <- list(list(X = GENO, model = 'BayesB'))
                  MODEL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = F)
                }
                
                gpred <- MODEL$yHat
                gpredSD <- MODEL$SD.yHat
                VARU <- var(gpred)
                SDU <- MODEL$SD.yHat
                
                CD <- sqrt(1 - (SDU^2 / VARU))
                GPRED <- cbind(data.frame(data.frame(gpred)[match(rownames(P_Test), rownames(data.frame(gpred))),]), data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BB Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "BC"){
                
                GENO <- rbind(G_Train, G_pred)
                colnames(P_Train) <- "V"
                colnames(P_Test) <- "V"
                y <- rbind(P_Train, P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa <- setNames(y$V, rownames(y))
                nIter <- 1500
                burnIn <- 500
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred <- as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain, FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2 <- list(list(X = FIX, model = "FIXED"), list(X = GENO, model = 'BayesC'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA2, nIter = nIter, burnIn = burnIn, verbose = F)
                } else {
                  options(warn = -1)
                  ETA <- list(list(X = GENO, model = 'BayesC'))
                  MODEL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = F)
                }
                
                gpred <- MODEL$yHat
                gpredSD <- MODEL$SD.yHat
                VARU <- var(gpred)
                SDU <- MODEL$SD.yHat
                
                CD <- sqrt(1 - (SDU^2 / VARU))
                GPRED <- cbind(data.frame(data.frame(gpred)[match(rownames(P_Test), rownames(data.frame(gpred))),]), data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BC Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "BL"){
                
                GENO <- rbind(G_Train, G_pred)
                colnames(P_Train) <- "V"
                colnames(P_Test) <- "V"
                y <- rbind(P_Train, P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa <- setNames(y$V, rownames(y))
                nIter <- 1500
                burnIn <- 500
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred <- as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain, FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2 <- list(list(X = FIX, model = "FIXED"), list(X = GENO, model = 'BL'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA2, nIter = nIter, burnIn = burnIn, verbose = F)
                } else {
                  options(warn = -1)
                  ETA <- list(list(X = GENO, model = 'BL'))
                  MODEL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = F)
                }
                
                gpred <- MODEL$yHat
                gpredSD <- MODEL$SD.yHat
                VARU <- var(gpred)
                SDU <- MODEL$SD.yHat
                
                CD <- sqrt(1 - (SDU^2 / VARU))
                GPRED <- cbind(data.frame(data.frame(gpred)[match(rownames(P_Test), rownames(data.frame(gpred))),]), data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BL Fold ', i, " Time ", n, ' finished #'))
              }
              
              #BayesRR
              if (pmethod == "BRR"){
                
                GENO <- rbind(G_Train, G_pred)
                colnames(P_Train) <- "V"
                colnames(P_Test) <- "V"
                y <- rbind(P_Train, P_Test)
                y[rownames(P_Test),] <- NA
                
                yNa <- setNames(y$V, rownames(y))
                nIter <- 1500
                burnIn <- 500
                
                if (Fixed_effect[1] != "NULL"){
                  FixedTrain <- as.matrix(FIXED[rownames(FIXED) %in% rownames(P_Train),])
                  FixedPred <- as.matrix(FIXED[!rownames(FIXED) %in% rownames(P_Train),])
                  
                  FIX <- data.frame(rbind(FixedTrain, FixedPred))
                  FIX <- FIX %>% mutate_if(is.character, as.numeric)
                  ETA2 <- list(list(X = FIX, model = "FIXED"), list(X = GENO, model = 'BRR'))
                  options(warn = -1)
                  MODEL <- BGLR(y = yNa, ETA = ETA2, nIter = nIter, burnIn = burnIn, verbose = F)
                } else {
                  options(warn = -1)
                  ETA <- list(list(X = GENO, model = 'BRR'))
                  MODEL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = F)
                }
                
                gpred <- MODEL$yHat
                gpredSD <- MODEL$SD.yHat
                VARU <- var(gpred)
                SDU <- MODEL$SD.yHat
                
                CD <- sqrt(1 - (SDU^2 / VARU))
                GPRED <- cbind(data.frame(data.frame(gpred)[match(rownames(P_Test), rownames(data.frame(gpred))),]), data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BRR Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "EN"){
                
                cv.fit <- cv.glmnet(G_Train, P_Train[,1], family = "gaussian", alpha = 0.5, nfolds = 10)
                lambda_min <- cv.fit$lambda.min
                gpred <- predict(cv.fit, newx = G_pred, s = c(lambda_min))
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('EN Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "BRNN"){
                
                cv.fit <- brnn(G_Train, P_Train[,1], neurons = 2, epochs = 10, verbose = T)
                gpred <- predict(cv.fit, newdata = G_pred)
                
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('BRNN Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "RR"){
                
                cv.fit <- cv.glmnet(G_Train, P_Train[,1], alpha = 0)
                lambda_min <- cv.fit$lambda.min
                gpred <- predict(cv.fit, newx = G_pred, s = c(lambda_min))
                
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('RR Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "RF"){
                
                gpred <- randomForest(G_Train, P_Train[,1], xtest = G_pred)$test$predicted
                
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('RF Fold ', i, " Time ", n, ' finished #'))
              }
              
              if (pmethod == "SVM"){
                
                model <- svm(G_Train, P_Train, method = "nu-regression", kernel = "radial", cost = 10, gamma = 0.001)
                gpred <- predict(model, G_pred)
                
                gpredSD <- NA
                CD <- NA
                GPRED <- cbind(gpred, data.frame(P_Test))
                GPRED <- round(GPRED, digits = 3)
                rownames(GPRED) <- rownames(G_pred)
                colnames(GPRED) <- c(paste0('Predicted', n), paste0('Actual', n))
                message(paste('SVM Fold ', i, " Time ", n, ' finished #'))
              }
              
              GPRED <- as.matrix(GPRED)
              df_total <- rbind(df_total, GPRED)
            }
            
            df[[n]] <- df_total
          }
          
          df_t <- data.frame()
          for (r in 1:iter){
            if (r == 1){
              df_t <- data.frame(df[[r]])
            } else {
              df_m <- data.frame(df[[r]])
              df_mt <- df_m[match(rownames(df_t), rownames(df_m)),]
              df_t <- cbind(df_t, df_mt)
            }
          }
          
          predicted_cols <- grep("Predicted", names(df_t))
          actual_cols <- grep("Actual", names(df_t))
          
          # Calculate and add mean columns
          df_t$PredictedMean <- rowMeans(df_t[, predicted_cols])
          df_t$ActualMean <- rowMeans(df_t[, actual_cols])
          
          # Round the mean columns
          df_t[, c("PredictedMean", "ActualMean")] <- round(df_t[, c("PredictedMean", "ActualMean")], digits = 3)
          df1 <- df_t
          write.table(df1, paste(data_path, pmethod, "_Folds_", i, "_Times_", n, "_", trait, "_", 
                                 paste(Fixed_effect, collapse = "_"), "_", popsize, "_", markersize, ".txt", sep = ""), 
                      quote = F, sep = "\t")
          print(cor(df1$PredictedMean, df1$ActualMean))
          accuracy_results[[paste(trait, pmethod, popsize, markersize)]] <- cor(df1$PredictedMean, df1$ActualMean)
          
        } #end of trait names
      } #end of pmethods
    } #end of marker sizes
  } #end of pop sizes
  
  Accuracy <- unlist(accuracy_results)
  write.table(Accuracy, paste0(data_path, "accuracy_results", ".txt"), sep = "\t", quote = F, col.names = F)
}

################################################################################
# EXAMPLE USAGE:  Compare regular GBLUP vs Weighted GBLUP
################################################################################

# 1.Run regular GBLUP
run_CV(pheno_file = "pheno",
       geno_file = "M",
       fold = 5,
       iter = 10,
       trait_names = c("Is_Dead", "DPC_date", "DPC_time"),
       pmethods = c("GBLUP"),
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       marker_sizes = c("top_1000", "top_5000", "top_10000", "ALL"),
       data_path = paste0("data", "/")
)

# 2.Run Weighted GBLUP with inverse p-value weighting
run_CV(pheno_file = "pheno",
       geno_file = "M",
       fold = 5,
       iter = 10,
       trait_names = c("DPC_time"),
       pmethods = c("WGBLUP"),  # NEW METHOD
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       marker_sizes = c("50", "100", "20000", "30000", "40000", "ALL"),
       data_path = paste0("data", "/"),
       weight_method = "inverse_pvalue",  # Options: "inverse_pvalue", "sqrt_inverse", "power", "bonferroni"
       weight_alpha = 1  # Only used if weight_method = "power"
)

# 3.Compare both methods side-by-side
run_CV(pheno_file = "pheno",
       geno_file = "M",
       fold = 5,
       iter = 10,
       trait_names = c("DPC_time"),
       pmethods = c("GBLUP", "WGBLUP"),  # Both methods
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       marker_sizes = c("top_1000"),
       data_path = paste0("data", "/"),
       weight_method = "inverse_pvalue",
       weight_alpha = 1
)

# 4.Test different weighting strategies
run_CV(pheno_file = "pheno",
       geno_file = "M",
       fold = 5,
       iter = 10,
       trait_names = c("DPC_date"),
       pmethods = c("WGBLUP"),
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       marker_sizes = c("top_5000"),
       data_path = paste0("data", "/"),
       weight_method = "sqrt_inverse"  # Try different weighting
)

# 
###############################################################################
# Plotting WGBLUP vs GBLUP Comparison with Friedman Test and Ranking Letters
################################################################################

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
    cat("\nFriedman test is significant (p < 0.05). Performing post-hoc analysis.. .\n")
    
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
      # Fallback:  simple ranking based on means
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
  
  cat("Generating ranking letters.. .\n")
  cat("Number of methods:", n_methods, "\n")
  
  # Check if p_matrix is valid
  if(is.null(p_matrix) || ! is.matrix(p_matrix)) {
    cat("Invalid p-value matrix.Using mean-based ranking.\n")
    return(assign_letters_by_mean(mean_ranks))
  }
  
  # Get matrix row and column names
  p_rownames <- rownames(p_matrix)
  p_colnames <- colnames(p_matrix)
  
  # Check if all methods are present in the matrix
  methods_in_matrix <- methods[methods %in% p_rownames & methods %in% p_colnames]
  
  if(length(methods_in_matrix) < length(methods)) {
    cat("Warning: Not all methods found in p-value matrix.\n")
    cat("Using mean-based ranking as fallback.\n")
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
          if(! sig_diff[method, other_method]) {
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
  
  # Simple ranking:  each method gets a different letter based on rank
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

################################################################################
# MAIN ANALYSIS:  Compare WGBLUP vs GBLUP across marker sizes
################################################################################

# Parse Arguments
dir <- "data"
trait_name <- "DPC_time"  # Change to your trait
iter <- 10
fold <- 5
pop_size <- "NULL"
Fixed_effect <- "NULL"
h2 <- 1

# Define methods to compare
methods <- c("GBLUP", "WGBLUP")
methods <- c()

# Define marker sizes to test
marker_sizes <- c("top_500","top_1000", "top_5000", "top_10000", "ALL")
# Or use: marker_sizes <- c("top_50", "top_100", "top_500", "top_1000", "top_5000", "top_10000", "top_20000", "top_30000", "top_40000", "top_50000", "ALL")

data_path <- paste0(dir, "/")

# Create comparison matrix for all combinations
total_configs <- length(marker_sizes) * length(methods)
compareM <- matrix(0, iter, total_configs)

# Create configuration names and track which are which
config_names <- c()
config_counter <- 1

# Process each combination
for (marker_idx in seq_along(marker_sizes)) {
  for (method_idx in seq_along(methods)) {
    marker_size <- marker_sizes[marker_idx]
    pmethod <- methods[method_idx]
    
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, 
                       "_", trait_name, "_", Fixed_effect, "_", pop_size, 
                       "_", marker_size, ".txt")
    
    cat("Looking for file:", filename, "\n")
    
    # Check if file exists
    if (file.exists(filename)) {
      cat("File found, reading data...\n")
      datap <- read.table(filename, header = TRUE)
      
      for (n in 1:iter) {
        pred_col <- paste0("Predicted", n)
        actual_col <- paste0("Actual", n)
        
        if (pred_col %in% colnames(datap) && actual_col %in% colnames(datap)) {
          compareM[n, config_counter] <- cor(datap[, pred_col], datap[, actual_col])
        }
      }
    } else {
      cat("Warning: File not found:", filename, "\n")
      compareM[, config_counter] <- NA
    }
    
    # Create configuration name
    config_names <- c(config_names, paste0(pmethod, "_", marker_size))
    config_counter <- config_counter + 1
  }
}

# Set column and row names
colnames(compareM) <- config_names
rownames(compareM) <- paste0("CV", 1:iter)

# Remove columns with all NA values
compareM <- compareM[, ! apply(is.na(compareM), 2, all)]

# Check if we have valid data
if (ncol(compareM) == 0) {
  stop("No valid data files found.Please check file paths and names.")
}

namef <- paste0("_", trait_name, "_WGBLUP_vs_GBLUP_markers")

# Save the comparison matrix
write.table(compareM, paste0(data_path, "Table", namef, ".csv"), 
            sep = ",", quote = FALSE)

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
write.csv(friedman_summary, paste0(data_path, "Friedman_Results", namef, ".csv"), 
          row.names = FALSE)

# Create ranking summary table
ranking_summary <- data.frame(
  Configuration = names(ranking_letters),
  Mean_Prediction_Ability = round(colMeans(compareM, na.rm = TRUE), 4),
  SD = round(apply(compareM, 2, sd, na.rm = TRUE), 4),
  Ranking_Letter = ranking_letters[names(ranking_letters)]
) %>%
  arrange(desc(Mean_Prediction_Ability))

write.csv(ranking_summary, paste0(data_path, "Rankings", namef, ".csv"), 
          row.names = FALSE)

################################################################################
# CREATE VISUALIZATION
################################################################################

# Melt the data for plotting
plot_data <- reshape2::melt(compareM)
colnames(plot_data) <- c("CV_Iteration", "Configuration", "Predictive_Ability")

# Extract method and marker size
plot_data$Method <- gsub("_top_.*|_ALL.*", "", plot_data$Configuration)
plot_data$Marker_Size <- gsub(".*_(top_[0-9]+|ALL).*", "\\1", plot_data$Configuration)

# Extract numeric values for ordering
plot_data$numeric_value <- gsub("top_", "", plot_data$Marker_Size)
plot_data$numeric_value[plot_data$numeric_value == "ALL"] <- "999999"  # Large number for sorting
plot_data$numeric_value <- as.numeric(plot_data$numeric_value)

# Create proper ordering for x-axis
marker_order_numeric <- sort(unique(plot_data$numeric_value))
marker_labels <- ifelse(marker_order_numeric == 999999, "ALL", 
                        as.character(marker_order_numeric))

plot_data$x_label <- factor(plot_data$numeric_value, 
                            levels = marker_order_numeric,
                            labels = marker_labels)

# Calculate y-positions for ranking letters
max_values <- plot_data %>%
  group_by(x_label, Method) %>%
  summarise(max_val = max(Predictive_Ability, na.rm = TRUE), .groups = 'drop')

# Add ranking letters to max_values
max_values$ranking_letter <- ""
for (i in 1:nrow(max_values)) {
  marker_label <- as.character(max_values$x_label[i])
  method_val <- max_values$Method[i]
  
  # Reconstruct configuration name
  if (marker_label == "ALL") {
    config_name <- paste0(method_val, "_ALL")
  } else {
    config_name <- paste0(method_val, "_top_", marker_label)
  }
  
  if (config_name %in% names(ranking_letters)) {
    max_values$ranking_letter[i] <- ranking_letters[config_name]
  }
}

# Define colors
method_colors <- c("GBLUP" = "#443A83FF", "WGBLUP" = "#27AD81FF")

# Create the first plot - Prediction Accuracy (with h2 correction)
p1 <- ggplot(plot_data, aes(x = x_label, y = Predictive_Ability/sqrt(h2), fill = Method)) + 
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = method_colors,
    name = "Method",
    labels = c("GBLUP", "Weighted GBLUP")
  ) +
  labs(title = paste0("Prediction Accuracy - ", trait_name), 
       x = "Number of Markers", 
       y = "Prediction Accuracy") +
  # Add mean points
  stat_summary(aes(group = Method), fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  # Add ranking letters
  geom_text(data = max_values, 
            aes(x = x_label, y = (max_val/sqrt(h2)) + 0.03, 
                label = ranking_letter, group = Method),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", color = "blue", inherit.aes = FALSE) +
  theme(
    legend.title = element_text(colour = "blue", size = 16, face = "bold", margin = margin(b = 15)),
    legend.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "#FFFFFF"),
    legend.background = element_rect(linetype = "dashed"),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(data_path, "Accuracy_Plot", namef, ".png"), 
       p1, scale = 1, dpi = 300, height = 6, width = 12)

# Create the second plot - Predictive Ability (without h2 correction)
p2 <- ggplot(plot_data, aes(x = x_label, y = Predictive_Ability, fill = Method)) + 
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = method_colors,
    name = "Method",
    labels = c("GBLUP", "Weighted GBLUP")
  ) +
  labs(title = paste0(trait_name), 
       x = "Number of Markers", 
       y = "Predictive Ability") +
  # Add mean points
  stat_summary(aes(group = Method), fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  # Add ranking letters
  geom_text(data = max_values, 
            aes(x = x_label, y = max_val + 0.03, 
                label = ranking_letter, group = Method),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", color = "blue", inherit.aes = FALSE) +
  theme(
    legend.title = element_text(colour = "blue", size = 16, face = "bold", margin = margin(b = 15)),
    legend.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "#FFFFFF"),
    legend.background = element_rect(linetype = "dashed"),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(data_path, "Ability_Plot", namef, ".png"), 
       p2, scale = 1, dpi = 300, height = 6, width = 12)

################################################################################
# PRINT SUMMARY
################################################################################

cat("\n=== WGBLUP vs GBLUP MARKER COMPARISON ANALYSIS COMPLETE ===\n")
cat("Files generated:\n")
cat("1.Table", namef, ".csv\n")
cat("2.Accuracy_Plot", namef, ".png\n")
cat("3.Ability_Plot", namef, ".png\n")
cat("4.Friedman_Results", namef, ". csv\n")
cat("5.Rankings", namef, ".csv\n")

cat("\n=== RANKING SUMMARY ===\n")
print(ranking_summary)

cat("\nConfiguration Ranking Letters:\n")
for(i in 1:length(ranking_letters)) {
  cat(names(ranking_letters)[i], ": ", ranking_letters[i], "\n")
}

cat("\nNote: Configurations sharing the same letter are not significantly different.\n")
cat("Letters are assigned alphabetically with 'a' being the best performing group.\n")

################################################################################
# Plotting GBLUP-GWAS vs WGBLUP-Random vs GBLUP-Random Comparison with Friedman Test and Ranking Letters
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
    
    if(! is.null(nemenyi_result)) {
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
      # Fallback:  simple ranking based on means
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
  
  cat("Generating ranking letters.. .\n")
  cat("Number of methods:", n_methods, "\n")
  
  # Check if p_matrix is valid
  if(is.null(p_matrix) || ! is.matrix(p_matrix)) {
    cat("Invalid p-value matrix.Using mean-based ranking.\n")
    return(assign_letters_by_mean(mean_ranks))
  }
  
  # Get matrix row and column names
  p_rownames <- rownames(p_matrix)
  p_colnames <- colnames(p_matrix)
  
  # Check if all methods are present in the matrix
  methods_in_matrix <- methods[methods %in% p_rownames & methods %in% p_colnames]
  
  if(length(methods_in_matrix) < length(methods)) {
    cat("Warning: Not all methods found in p-value matrix.\n")
    cat("Using mean-based ranking as fallback.\n")
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
          if(! sig_diff[method, other_method]) {
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
  
  # Simple ranking:  each method gets a different letter based on rank
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

################################################################################
# MAIN ANALYSIS:  Compare GBLUP-GWAS vs WGBLUP-Random vs GBLUP-Random
################################################################################

# Parse Arguments
dir <- "data"
trait_name <- "DPC_time"  # Change to your trait
iter <- 10
fold <- 5
pop_size <- "NULL"
Fixed_effect <- "NULL"
h2 <- 1

# Define marker sizes to test
# For GWAS markers (GBLUP): "top_X"
# For Random markers (GBLUP and WGBLUP): numeric "X"
marker_sizes_gwas <- c("top_50", "top_100", "top_1000", "top_5000", "top_10000", "top_20000", "top_30000", "top_40000", "top_50000")
marker_sizes_random <- c("50", "100", "1000", "5000", "10000", "20000", "30000", "40000", "50000")

# Extract numeric values for x-axis labels
marker_numbers <- c("1000", "5000", "10000", "20000", "30000", "40000", "50000")

data_path <- paste0(dir, "/")

# Create comparison matrix for all combinations
total_configs <- length(marker_numbers) * 3  # 3 methods
compareM <- matrix(0, iter, total_configs)

# Create configuration names and track which are which
config_names <- c()
config_counter <- 1

# Process each marker size
for (marker_idx in seq_along(marker_numbers)) {
  marker_num <- marker_numbers[marker_idx]
  
  # 1.GBLUP with GWAS markers (top_X) - CYAN/GREEN
  pmethod <- "GBLUP"
  marker_size <- marker_sizes_gwas[marker_idx]
  
  filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, 
                     "_", trait_name, "_", Fixed_effect, "_", pop_size, 
                     "_", marker_size, ". txt")
  
  cat("Looking for GBLUP-GWAS file:", filename, "\n")
  
  if (file.exists(filename)) {
    cat("File found, reading data...\n")
    datap <- read.table(filename, header = TRUE)
    
    for (n in 1:iter) {
      pred_col <- paste0("Predicted", n)
      actual_col <- paste0("Actual", n)
      
      if (pred_col %in% colnames(datap) && actual_col %in% colnames(datap)) {
        compareM[n, config_counter] <- cor(datap[, pred_col], datap[, actual_col])
      }
    }
  } else {
    cat("Warning: File not found:", filename, "\n")
    compareM[, config_counter] <- NA
  }
  
  config_names <- c(config_names, paste0("GWAS_", marker_num))
  config_counter <- config_counter + 1
  
  # 2.WGBLUP with Random markers (numeric X) - ORANGE
  pmethod <- "WGBLUP"
  marker_size <- marker_sizes_random[marker_idx]
  
  filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, 
                     "_", trait_name, "_", Fixed_effect, "_", pop_size, 
                     "_", marker_size, ". txt")
  
  cat("Looking for WGBLUP-Random file:", filename, "\n")
  
  if (file.exists(filename)) {
    cat("File found, reading data...\n")
    datap <- read.table(filename, header = TRUE)
    
    for (n in 1:iter) {
      pred_col <- paste0("Predicted", n)
      actual_col <- paste0("Actual", n)
      
      if (pred_col %in% colnames(datap) && actual_col %in% colnames(datap)) {
        compareM[n, config_counter] <- cor(datap[, pred_col], datap[, actual_col])
      }
    }
  } else {
    cat("Warning: File not found:", filename, "\n")
    compareM[, config_counter] <- NA
  }
  
  config_names <- c(config_names, paste0("WGBLUPRandom_", marker_num))
  config_counter <- config_counter + 1
  
  # 3.GBLUP with Random markers (numeric X) - PURPLE
  pmethod <- "GBLUP"
  marker_size <- marker_sizes_random[marker_idx]
  
  filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, 
                     "_", trait_name, "_", Fixed_effect, "_", pop_size, 
                     "_", marker_size, ".txt")
  
  cat("Looking for GBLUP-Random file:", filename, "\n")
  
  if (file.exists(filename)) {
    cat("File found, reading data...\n")
    datap <- read.table(filename, header = TRUE)
    
    for (n in 1:iter) {
      pred_col <- paste0("Predicted", n)
      actual_col <- paste0("Actual", n)
      
      if (pred_col %in% colnames(datap) && actual_col %in% colnames(datap)) {
        compareM[n, config_counter] <- cor(datap[, pred_col], datap[, actual_col])
      }
    }
  } else {
    cat("Warning: File not found:", filename, "\n")
    compareM[, config_counter] <- NA
  }
  
  config_names <- c(config_names, paste0("Random_", marker_num))
  config_counter <- config_counter + 1
}

# Set column and row names
colnames(compareM) <- config_names
rownames(compareM) <- paste0("CV", 1:iter)

# Remove columns with all NA values
compareM <- compareM[, ! apply(is.na(compareM), 2, all)]

# Check if we have valid data
if (ncol(compareM) == 0) {
  stop("No valid data files found.Please check file paths and names.")
}

cat("\n=== COMPARISON MATRIX ===\n")
print(head(compareM))
cat("\n")

namef <- paste0("_", trait_name, "_GWAS_vs_WGBLUPRandom_vs_Random")

# Save the comparison matrix
write.table(compareM, paste0(data_path, "Table", namef, ".csv"), 
            sep = ",", quote = FALSE)

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
write.csv(friedman_summary, paste0(data_path, "Friedman_Results", namef, ".csv"), 
          row.names = FALSE)

# Create ranking summary table
ranking_summary <- data.frame(
  Configuration = names(ranking_letters),
  Mean_Prediction_Ability = round(colMeans(compareM, na.rm = TRUE), 4),
  SD = round(apply(compareM, 2, sd, na.rm = TRUE), 4),
  Ranking_Letter = ranking_letters[names(ranking_letters)]
) %>%
  arrange(desc(Mean_Prediction_Ability))

write.csv(ranking_summary, paste0(data_path, "Rankings", namef, ".csv"), 
          row.names = FALSE)

################################################################################
# CREATE VISUALIZATION
################################################################################

# Melt the data for plotting
plot_data <- reshape2::melt(compareM)
colnames(plot_data) <- c("CV_Iteration", "Configuration", "Predictive_Ability")

# Extract method and marker size
plot_data$Method <- gsub("_[0-9]+$", "", plot_data$Configuration)
plot_data$Marker_Size <- as.numeric(gsub(".*_([0-9]+)$", "\\1", plot_data$Configuration))

# Create proper ordering for x-axis
plot_data$x_label <- factor(plot_data$Marker_Size, 
                            levels = sort(unique(plot_data$Marker_Size)))

# Calculate y-positions for ranking letters
max_values <- plot_data %>%
  group_by(x_label, Method) %>%
  summarise(max_val = max(Predictive_Ability, na.rm = TRUE), .groups = 'drop')

# Add ranking letters to max_values
max_values$ranking_letter <- ""
for (i in 1:nrow(max_values)) {
  marker_label <- as.character(max_values$x_label[i])
  method_val <- max_values$Method[i]
  
  # Reconstruct configuration name
  config_name <- paste0(method_val, "_", marker_label)
  
  if (config_name %in% names(ranking_letters)) {
    max_values$ranking_letter[i] <- ranking_letters[config_name]
  }
}

# Define colors matching your specification
method_colors <- c(
  "GWAS" = "#27AD81FF",            # Cyan/Green (GBLUP with GWAS markers)
  "WGBLUPRandom" = "#FFA500",      # Orange (WGBLUP with Random markers)
  "Random" = "#443A83FF"           # Purple (GBLUP with Random markers)
)

# Create the first plot - Prediction Accuracy (with h2 correction)
p1 <- ggplot(plot_data, aes(x = x_label, y = Predictive_Ability/sqrt(h2), fill = Method)) + 
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = method_colors,
    name = "Marker Type",
    labels = c("GWAS" = "GBLUP-GWAS", "WGBLUPRandom" = "WGBLUP-Random", "Random" = "GBLUP-Random")
  ) +
  labs(title = paste0("Prediction Accuracy - ", trait_name), 
       x = "Number of Markers", 
       y = "Prediction Accuracy") +
  # Add mean points
  stat_summary(aes(group = Method), fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  # Add ranking letters
  geom_text(data = max_values, 
            aes(x = x_label, y = (max_val/sqrt(h2)) + 0.03, 
                label = ranking_letter, group = Method),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", color = "blue", inherit.aes = FALSE) +
  theme(
    legend.title = element_text(colour = "blue", size = 16, face = "bold", margin = margin(b = 15)),
    legend.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "#FFFFFF"),
    legend.background = element_rect(linetype = "dashed"),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(data_path, "Accuracy_Plot", namef, ".png"), 
       p1, scale = 1, dpi = 300, height = 6, width = 12)

# Create the second plot - Predictive Ability (without h2 correction)
p2 <- ggplot(plot_data, aes(x = x_label, y = Predictive_Ability, fill = Method)) + 
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = method_colors,
    name = "Marker Type",
    labels = c("GWAS" = "GBLUP-GWAS", "WGBLUPRandom" = "WGBLUP-Random", "Random" = "GBLUP-Random")
  ) +
  labs(title = paste0(trait_name), 
       x = "Number of Markers", 
       y = "Predictive Ability") +
  # Add mean points
  stat_summary(aes(group = Method), fun = mean, geom = "point", shape = 21, size = 2, 
               color = "black", fill = "red", stroke = 1,
               position = position_dodge(width = 0.8)) +
  # Add ranking letters
  geom_text(data = max_values, 
            aes(x = x_label, y = max_val + 0.03, 
                label = ranking_letter, group = Method),
            position = position_dodge(width = 0.8),
            size = 5, fontface = "bold", color = "blue", inherit.aes = FALSE) +
  theme(
    legend.title = element_text(colour = "blue", size = 16, face = "bold", margin = margin(b = 15)),
    legend.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    plot.background = element_rect(fill = "#FFFFFF"),
    legend.background = element_rect(linetype = "dashed"),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(data_path, "Ability_Plot", namef, ".png"), 
       p2, scale = 1, dpi = 300, height = 6, width = 12)

################################################################################
# PRINT SUMMARY
################################################################################

cat("\n=== GBLUP-GWAS vs WGBLUP-Random vs GBLUP-Random COMPARISON COMPLETE ===\n")
cat("Files generated:\n")
cat("1.Table", namef, ".csv\n")
cat("2.Accuracy_Plot", namef, ".png\n")
cat("3.Ability_Plot", namef, ". png\n")
cat("4.Friedman_Results", namef, ".csv\n")
cat("5.Rankings", namef, ".csv\n")

cat("\n=== RANKING SUMMARY ===\n")
print(ranking_summary)

cat("\nConfiguration Ranking Letters:\n")
for(i in 1:length(ranking_letters)) {
  cat(names(ranking_letters)[i], ": ", ranking_letters[i], "\n")
}

cat("\nNote:  Configurations sharing the same letter are not significantly different.\n")
cat("Letters are assigned alphabetically with 'a' being the best performing group.\n")

# Additional pairwise comparisons between methods
cat("\n=== PAIRWISE METHOD COMPARISONS ===\n")
for(marker_num in marker_numbers) {
  cat("\n--- Marker Size:", marker_num, "---\n")
  
  # Get data for this marker size
  gwas_col <- paste0("GWAS_", marker_num)
  wgblup_random_col <- paste0("WGBLUPRandom_", marker_num)
  random_col <- paste0("Random_", marker_num)
  
  if(all(c(gwas_col, wgblup_random_col, random_col) %in% colnames(compareM))) {
    gwas_data <- compareM[, gwas_col]
    wgblup_random_data <- compareM[, wgblup_random_col]
    random_data <- compareM[, random_col]
    
    # GBLUP-GWAS vs WGBLUP-Random
    if(! any(is.na(gwas_data)) && !any(is.na(wgblup_random_data))) {
      t_test <- t.test(gwas_data, wgblup_random_data, paired = TRUE)
      cat("GBLUP-GWAS vs WGBLUP-Random:  p-value =", round(t_test$p.value, 4), 
          ifelse(t_test$p.value < 0.05, " *SIGNIFICANT*", ""), "\n")
      cat("  Mean difference:", round(mean(gwas_data) - mean(wgblup_random_data), 4), "\n")
    }
    
    # GBLUP-GWAS vs GBLUP-Random
    if(! any(is.na(gwas_data)) && !any(is.na(random_data))) {
      t_test <- t.test(gwas_data, random_data, paired = TRUE)
      cat("GBLUP-GWAS vs GBLUP-Random: p-value =", round(t_test$p.value, 4), 
          ifelse(t_test$p.value < 0.05, " *SIGNIFICANT*", ""), "\n")
      cat("  Mean difference:", round(mean(gwas_data) - mean(random_data), 4), "\n")
    }
    
    # WGBLUP-Random vs GBLUP-Random
    if(!any(is.na(wgblup_random_data)) && !any(is.na(random_data))) {
      t_test <- t.test(wgblup_random_data, random_data, paired = TRUE)
      cat("WGBLUP-Random vs GBLUP-Random: p-value =", round(t_test$p.value, 4), 
          ifelse(t_test$p.value < 0.05, " *SIGNIFICANT*", ""), "\n")
      cat("  Mean difference:", round(mean(wgblup_random_data) - mean(random_data), 4), "\n")
    }
  }
}

