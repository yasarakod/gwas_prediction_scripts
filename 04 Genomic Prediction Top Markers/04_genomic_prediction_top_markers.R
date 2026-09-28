# ============================================================================== #
# Script: 04_genomic_prediction_top_markers.R
# Title: Genomic Prediction with Top GWAS-Selected Markers
# Description:
#   Executes 10-fold cross-validation for target traits (DPC_time, DPC_date, Is_Dead) utilizing prioritized top GWAS marker sets across all Bayesian, BLUP, and machine learning models.
#
# Inputs:
#   - data/M (Genotype matrix)
#   - data/pheno (Phenotype records)
#   - data/top_1000_snps.txt (Top 1000 GWAS markers selected from Wald test)
#
# Outputs:
#   - results/*_Folds_5_Times_10_*_top_1000.txt (Cross-validation prediction logs for top 1000 markers)
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
          
          # Assuming your file is named "snp_pvalues.txt" and has two columns: "snp_id" and "p_value"
          snp_data <- read.table(paste0(data_path, "wald_",trait,".txt"), header = TRUE)
          
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
          
          # Function to randomly sample markers
          # geno_shrink <- ped_formated[, if (!is.null(markersize) && markersize < ncol(ped_formated)) 
          #   sample(ncol(ped_formated), markersize) 
          #   else 1:ncol(ped_formated), drop = FALSE]
          
          # Function to selectively sample markers or use top SNPs from files
          
          geno_data_shrink <- function(ped_formated, markersize, data_path) {
            if (is.character(markersize) && grepl("^top_", markersize)) {
              # Read top SNPs from file
              top_snps <- read.table(paste0(data_path, markersize, "_snps.txt"), header = FALSE)$V1
              top_snps <- gsub("-", ".", top_snps)
              
              # Subset based on top SNPs
              geno_data_shrink <- ped_formated[, colnames(ped_formated) %in% top_snps, drop = FALSE]
              
              message("Using top SNPs for marker size: ", markersize)
            } else if (markersize == "ALL") {
              # Use all markers
              geno_data_shrink <- ped_formated
              
              message("Using all markers.")
            } else {
              # Randomly sample markers if markersize is numeric
              geno_data_shrink <- ped_formated[, if (!is.null(markersize) && markersize < ncol(ped_formated)) 
                sample(ncol(ped_formated), markersize) 
                else 1:ncol(ped_formated), drop = FALSE] 
            }
            
            return(geno_data_shrink)
          }
          # Call the geno_shrink function and assign the output to a new variable
          geno_shrink <- geno_data_shrink(ped_formated, markersize, data_path)
          
          # Use the new variable to print the dimensions
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
              
              #GBLUP
              if (pmethod == "GBLUP"){
                
                if(Fixed_effect[1]!="NULL"){
                  #saveRDS(FIXED,"FIXED")
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
  write.table(Accuracy,paste0(data_path,"accuracy_results",".txt"),sep = "\t",quote = F, col.names = F)
}

run_CV(pheno_file="pheno",
       geno_file="M",
       fold = 5,
       iter = 10,
       trait_names = c("Is_Dead","DPC_date","DPC_time"),
       pmethods= c("PBLUP","GBLUP","BA","BB","BC","BRR","BRNN","RR","EN","RF"),
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       marker_sizes = c("NULL"),
       #pop_sizes = c(50, 100, 150, 200, 250, 300, 350, 400, 474),
       # marker_sizes = c("top_50","top_60","top_70","top_80","top_90","top_100",
       #                  "top_500","top_1000","top_5000","top_10000","top_20000",
       #                  "top_30000","top_40000","top_50000","ALL"),
       #marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL"),
       data_path = paste0("data", "/")
)

run_CV(pheno_file="pheno",
       geno_file="M",
       fold = 5,
       iter = 10,
       trait_names = c("Is_Dead","DPC_date","DPC_time"),
       #trait_names = c("Is_Dead","DPC_date","DPC_time"),
       pmethods= c("GBLUP"),
       Fixed_effect = c("NULL"),
       pop_sizes = c("NULL"),
       # marker_sizes = c("top_50","top_60","top_70","top_80","top_90","top_100",
       #                  "top_500","top_1000","top_5000","top_10000","top_20000",
       #                  "top_30000","top_40000","top_50000","ALL"),
       marker_sizes = c("top_1000"),
       data_path = paste0("data", "/")
)

run_CV(pheno_file="pheno",
       geno_file="M",
       fold = 5,
       iter = 10,
       trait_names = c("Is_Dead","DPC_date","DPC_time"),
       pmethods= c("GBLUP"),
       Fixed_effect = c("NULL"),
       pop_sizes = c("50","100", "150", "200", "250", "300", "350", "400", "474"),
       marker_sizes = c("top_1000"),
       data_path = paste0("data", "/")
)

import_packages<-function(){
  #define packages to install
  packages <- c("argparse","ggplot2","reshape2","ggpubr")
  
  install.packages(setdiff(packages, rownames(installed.packages())))
  
  invisible(lapply(packages, library, character.only=TRUE))
}

import_packages()

# parser <- ArgumentParser(description = "CV graphs")
# 
# # Define Arguments 
# parser$add_argument("--dir", default = "data", required = TRUE,help = "Data directory")
# parser$add_argument("--graph_type", required = TRUE, help = "Graph type")
# parser$add_argument("--trait_name",required = TRUE, help = "Trait name")
# parser$add_argument("--pmethods", nargs = "+", default = c("GBLUP","BA","BB","BC","BL","BRR","EN","RR","RF"), help = "Prediction methods (space-separated)")
# parser$add_argument("--iter", type = "integer", default = 3, help = "Number of iterations")
# parser$add_argument("--fold", type = "integer", default = 3, help = "Number of folds")
# parser$add_argument("--pop_sizes", nargs = "+", default = c("NULL"), help = "Population sizes (space-separated)")
# parser$add_argument("--marker_sizes", nargs = "+", default = c("NULL"), help = "Marker sizes (space-separated)")
# parser$add_argument("--Fixed_effect", nargs = "+", default = c("NULL"), help = "Fixed effects (space-separated)")
# parser$add_argument("--h2", type = "double",required = TRUE, help = "Heritability")
# #parser$add_argument("--marker_start", type = "integer", help = "Starting marker size")
# #parser$add_argument("--marker_end", type = "integer", help = "Ending marker size")


# Parse Arguments
dir<-"data"
graph_type<-"A"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("DPC_time")
pmethods<-c("PBLUP","GBLUP","BA","BB","BC","BL","BRR","RR","EN","RF")
iter<-10
fold<-5
pop_sizes<-"NULL"
#marker_sizes<-"NULL"
marker_sizes<-c("top_40000")
#marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("NULL")
h2<-1

data_path = paste0(dir, "/")

# Parse Arguments
dir<-"data"
graph_type<-"B"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("EP_Antibody_log")
pmethods<-c("BL")
iter<-10
fold<-5
pop_sizes<-"NULL"
#marker_sizes<-"NULL"
#marker_sizes<-c("top_40000")
marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("Plate")
h2<-1

data_path = paste0(dir, "/")

# Parse Arguments
dir<-"data"
graph_type<-"C"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("DPC_time")
pmethods<-c("GBLUP")
iter<-10
fold<-5
pop_sizes<-"NULL"
#marker_sizes<-"NULL"
marker_sizes<-c("top_1000")
#marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("NULL")
h2<-1

data_path = paste0(dir, "/")

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
  data<-reshape2::melt(compareM)
  
  table_graphic <- ggtexttable(round(compareM,4))
  ggsave(paste0(data_path, "Table", namef, ".png"), table_graphic, scale=1, dpi = 300, height = 2, width = 8)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.3)+
    labs(title=paste0("Plot_",trait),x= "Prediction method",fill="Method", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Plot",namef,".png"),p,scale=1, dpi = 300, height = 4, width = 4)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.3)+
    labs(title=paste0("Plot_",trait),x= "Prediction method",fill="Method", y = "Predictive Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Ability_Plot",namef,".png"),p,scale=1, dpi = 300, height = 4, width = 4)
}

if (graph_type=="B"){
  
  pmethod<-pmethods
  sname<-trait_name
  marker_sizes<-marker_sizes
  # Use the marker_sizes generated earlier
  compareM <- matrix(0, iter, length(marker_sizes))
  
  for (i in seq_along(marker_sizes)) {  
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, "_", sname, "_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_", marker_sizes[i], ".txt")
    datap <- read.table(filename)
    for (n in 1:iter) {
      compareM[n, i] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  # Use marker_sizes for column names
  colnames(compareM) <- paste0("m=", marker_sizes)
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_",sname,"_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_marker")
  
  write.table(compareM,paste0(data_path,"Table_",pmethod,namef,".csv"),sep = ",",quote = F)
  
  data<-reshape2::melt(compareM)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.37, 0.81)+
    labs(title=paste0("Plot_",sname),x= "Marker No",fill="Markers", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.37, 0.81)+
    labs(title=paste0("Plot_",sname),x= "Marker No",fill="Markers", y = "Prediction Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_text(angle = 45, size = 8),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Ability_Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
}
if (graph_type=="C"){
  
  sname<-trait_name
  pop_sizes<-pop_sizes
  pmethod<-pmethods
  compareM<-matrix(0,iter,length(pop_sizes))
  
  for (i in seq_along(pop_sizes)) {  
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, "_", sname, "_", Fixed_effect, "_", pop_sizes[i], "_", marker_sizes, ".txt") 
    datap <- read.table(filename)
    for (n in 1:iter) {
      compareM[n, i] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  colnames(compareM) <- paste0("s=", pop_sizes)
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_",sname,"_", Fixed_effect, "_pop_", marker_sizes)
  
  write.table(compareM,paste0(data_path,"Table_",pmethod,namef,".csv"),sep = ",",quote = F)
  
  data<-reshape2::melt(compareM)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.3, 0.8)+
    labs(title=paste0("Plot_",sname),x= "Pop size",fill="Pop size", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.3, 0.8)+
    labs(title=paste0("Plot_",sname),x= "Pop size",fill="Pop size", y = "Predictive Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Ability_Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
}



