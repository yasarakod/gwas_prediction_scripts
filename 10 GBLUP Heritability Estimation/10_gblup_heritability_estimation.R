# ============================================================================== #
# Script: 10_gblup_heritability_estimation.R
# Title: GBLUP Genomic Heritability Estimation
# Description:
#   Estimates genomic narrow-sense heritability (h2) for mortality and biometric traits using genomic relationship matrices via mixed model variance component estimation.
#
# Inputs:
#   - data/M (Genotype matrix)
#   - data/pheno (Phenotype observations: Is_Dead, DPC_time, DPC_date, Length, Weight, Width)
#
# Outputs:
#   - results/heritability_<trait>.rds (Estimated variance components and heritability estimates)
#   - results/heritability_summary.csv (Summary table of additive and residual variance)
#
# Required R Packages: rrBLUP, data.table, SNPRelate, snpReady
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


# Calculate Heritability from Your Processed Data Files
# Set working directory to your data path
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
  
  # # Read and filter the .fam file
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
  #detach("package:dplyr", unload = TRUE)
  base::system(paste0("plink.exe --allow-extra-chr --keep ",selection," --bfile ",infile," --make-bed --out ",data_path,"tempgeno"))
  base::system(paste0("plink.exe --allow-extra-chr --bfile ",data_path,"tempgeno --recode A --out ",data_path,"tempgeno"))
  base::system(paste0("plink.exe --allow-extra-chr --bfile ",data_path,"tempgeno --recode -tab --out ",data_path,"tempgeno"))
  
  ped<-data.frame(fread(paste0(data_path,"tempgeno.raw"),header = T)[,-c(1,3:6)])
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
    stop("Phenotype file not found: ", pheno_file)
  }
  
  pheno <- read.csv(paste0(data_path,pheno_file), header = TRUE)
  
  # Check for matching sample IDs
  if (!all(rownames(M) %in% pheno$Sample_ID)) {
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

#heritability calculation by rrBLUP

# Load required libraries
library(rrBLUP)

# Method 1: Using the final processed RDS files (RECOMMENDED)
calculate_heritability_from_rds <- function(trait_column = "trait_name") {
  
  # Load processed data
  cat("Loading processed data files...\n")
  M <- readRDS("data/M")          # Clean genotype matrix
  pheno <- readRDS("data/pheno")  # Matched phenotype data
  
  cat("Genotype matrix dimensions:", dim(M), "\n")
  cat("Phenotype data dimensions:", dim(pheno), "\n")
  cat("Available traits:", colnames(pheno), "\n")
  
  # Select trait for analysis
  if(trait_column %in% colnames(pheno)) {
    y <- pheno[, trait_column]
    cat("Analyzing trait:", trait_column, "\n")
  } else {
    cat("Available traits:", paste(colnames(pheno), collapse = ", "), "\n")
    stop("Trait '", trait_column, "' not found in phenotype data")
  }
  
  # Remove individuals with missing phenotypes
  complete_cases <- !is.na(y)
  M_clean <- M[complete_cases, ]
  y_clean <- y[complete_cases]
  
  cat("Complete observations:", sum(complete_cases), "\n")
  cat("Trait range:", round(range(y_clean), 3), "\n")
  
  # Create genomic relationship matrix
  cat("Creating genomic relationship matrix...\n")
  G <- A.mat(M_clean - 1)  # Convert 0,1,2 to -1,0,1 coding
  
  # Fit mixed model
  cat("Fitting GBLUP model...\n")
  model <- mixed.solve(y = y_clean, K = G)
  
  # Calculate heritability
  Vg <- model$Vu  # Genetic variance
  Ve <- model$Ve  # Error variance
  Vp <- Vg + Ve   # Phenotypic variance
  h2 <- Vg / Vp   # Narrow-sense heritability
  
  # Standard error
  n_obs <- length(y_clean)
  se_h2 <- sqrt((4 * h2 * (1 - h2)^2) / n_obs)
  
  # Print results
  cat("\n=== HERITABILITY RESULTS ===\n")
  cat("Trait:", trait_column, "\n")
  cat("Sample size:", n_obs, "\n")
  cat("Number of markers:", ncol(M_clean), "\n")
  cat("h² =", round(h2, 3), "± SE =", round(se_h2, 3), "\n")
  cat("Genetic variance (Vg) =", sprintf("%.6f", Vg), "\n")
  cat("Error variance (Ve) =", sprintf("%.6f", Ve), "\n")
  cat("Phenotypic variance (Vp) =", sprintf("%.6f", Vp), "\n")
  
  # Save results
  results <- list(
    trait = trait_column,
    h2 = h2,
    se_h2 = se_h2,
    Vg = Vg,
    Ve = Ve,
    Vp = Vp,
    n_obs = n_obs,
    n_markers = ncol(M_clean)
  )
  
  saveRDS(results, paste0("data/heritability_", trait_column, ".rds"))
  cat("Results saved to: data/heritability_", trait_column, ".rds\n")
  
  return(results)
}

# Method 2: Using original CSV files (if RDS files not available)
calculate_heritability_from_csv <- function(trait_column = "trait_name") {
  
  # Load original phenotype file
  pheno <- read.csv("data/phenotypes.csv", header = TRUE)
  cat("Original phenotype dimensions:", dim(pheno), "\n")
  cat("Available traits:", colnames(pheno), "\n")
  
  # Load selection file
  selected_samples <- read.table("data/selected.txt", header = FALSE)$V1
  cat("Selected samples:", length(selected_samples), "\n")
  
  # Filter phenotypes to selected samples
  pheno_filtered <- pheno[pheno$Sample_ID %in% selected_samples, ]
  cat("Filtered phenotype dimensions:", dim(pheno_filtered), "\n")
  
  # You would need to process genotypes here using your pipeline
  # This is more complex, so Method 1 is recommended
  
  cat("For this method, run your full preprocessing pipeline first,\n")
  cat("then use Method 1 with the resulting RDS files.\n")
}

# Method 3: Multiple traits analysis
analyze_all_traits <- function() {
  
  # Load data
  pheno <- readRDS("data/pheno")
  
  # Get all numeric columns (potential traits)
  trait_columns <- names(pheno)[sapply(pheno, is.numeric)]
  trait_columns <- trait_columns[trait_columns != "Sample_ID"]  # Exclude ID column
  
  cat("Found", length(trait_columns), "numeric traits:", paste(trait_columns, collapse = ", "), "\n\n")
  
  # Calculate heritability for each trait
  all_results <- list()
  
  for(trait in trait_columns) {
    cat("=== Analyzing trait:", trait, "===\n")
    tryCatch({
      result <- calculate_heritability_from_rds(trait)
      all_results[[trait]] <- result
    }, error = function(e) {
      cat("Error analyzing", trait, ":", e$message, "\n")
      all_results[[trait]] <- NULL
    })
    cat("\n")
  }
  
  # Summary table
  summary_df <- data.frame(
    Trait = names(all_results),
    h2 = sapply(all_results, function(x) if(!is.null(x)) round(x$h2, 3) else NA),
    SE = sapply(all_results, function(x) if(!is.null(x)) round(x$se_h2, 3) else NA),
    N_obs = sapply(all_results, function(x) if(!is.null(x)) x$n_obs else NA),
    stringsAsFactors = FALSE
  )
  
  cat("=== SUMMARY OF ALL TRAITS ===\n")
  print(summary_df)
  
  # Save summary
  write.csv(summary_df, "data/heritability_summary.csv", row.names = FALSE)
  saveRDS(all_results, "data/all_heritability_results.rds")
  
  return(all_results)
}

# Example usage:
cat("=== USAGE EXAMPLES ===\n")
cat("1.Single trait analysis:\n")
cat("   results <- calculate_heritability_from_rds('your_trait_name')\n\n")

cat("2.All traits analysis:\n") 
cat("   all_results <- analyze_all_traits()\n\n")

cat("3.Check available traits first:\n")
cat("   pheno <- readRDS('data/pheno')\n")
cat("   colnames(pheno)\n\n")

# Quick check function
check_data_availability <- function() {
  cat("=== DATA AVAILABILITY CHECK ===\n")
  
  files_to_check <- c(
    "data/M.rds",
    "data/pheno.rds", 
    "data/phenotypes.csv",
    "data/selected.txt",
    "data/hapmap.rds"
  )
  
  for(file in files_to_check) {
    if(file.exists(file)) {
      cat("✓", file, "- Available\n")
    } else {
      cat("✗", file, "- Missing\n")
    }
  }
  
  cat("\nTo generate missing RDS files, run your preprocessing pipeline:\n")
  cat("geno_data <- preprocess_genotype(...)\n")
  cat("geno_data <- perform_marker_qc(...)\n") 
  cat("geno_data <- prune_markers_by_ld(...)\n")
  cat("pheno_matched <- match_phenotype(...)\n")
}

# Run data check
check_data_availability()

#See available traits
pheno <- readRDS("data/pheno")
colnames(pheno)

# Calculate heritability for a specific trait
results <- calculate_heritability_from_rds("DPC_time")

# Or analyze all traits at once
all_results <- analyze_all_traits()
# Save all results
write.csv(all_results, "data/all_heritability_results.csv", row.names = FALSE)

