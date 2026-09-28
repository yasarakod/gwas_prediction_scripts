# Multi-Model Genomic Prediction Cross-Validation

**Script File:** `03_genomic_prediction_cross_validation.R`  
**Location:** `03 Genomic Prediction Cross Validation/03_genomic_prediction_cross_validation.R`

---

## 1. Brief Description
Executes k-fold cross-validation across 10 statistical, Bayesian, and machine learning prediction models (GBLUP, PBLUP, BayesA, BayesB, BayesC, BRR, BL, Elastic Net, Random Forest, rrBLUP) across various marker density subsets.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Required to generate the numerical genotype dosage matrix M, the formatted phenotype file pheno, and marker subsets top_*_snps.txt.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Required to compute the numerator additive pedigree relationship matrix A_matrix.rds needed for PBLUP benchmarking.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Preprocessed numeric genotype dosage matrix containing additive SNP alleles.
- **Expected File Format:** RDS serialized matrix or raw data matrix.
- **Header & Structure:** N rows (individuals) x M columns (SNP markers). Row names = individual IDs, Column names = marker IDs.
- **Value Types & Encodings:** Numeric integer or floating dosages in {0, 1, 2}. Missing values imputed.
- **How the Code Processes It:** Subsets to selected marker counts (e.g. top 50, 100, 500, 1000, ALL) and used as the predictor feature matrix in genomic prediction algorithms.

### 2. `pheno / pheno.rds`
- **Data Description:** Cleaned phenotype dataframe containing response variables.
- **Expected File Format:** Binary RDS file or dataframe.
- **Header & Structure:** Rows match exactly with rows of genotype matrix M. Columns include DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Continuous numeric values for time/date; 0/1 binary integers for survival.
- **How the Code Processes It:** Supplied as the target response vector Y in cross-validation training/validation splits.

### 3. `A_matrix.rds`
- **Data Description:** Numerator additive pedigree relationship matrix.
- **Expected File Format:** Binary RDS file containing N x N symmetric matrix.
- **Header & Structure:** Square covariance matrix with row/column names matching individual IDs.
- **Value Types & Encodings:** Floating-point additive kinship coefficients.
- **How the Code Processes It:** Used in pedigree BLUP (PBLUP) and combined single-step models.

### 4. `top_*_snps.txt (e.g. top_1000_snps.txt, top_500_snps.txt)`
- **Data Description:** Lists of prioritized SNP markers identified from GWAS Wald test significance rankings.
- **Expected File Format:** Single-column plain text file.
- **Header & Structure:** One marker identifier (SNP name) per line without header.
- **Value Types & Encodings:** Character strings matching column names of genotype matrix M.
- **How the Code Processes It:** Filters matrix M to evaluate genomic prediction accuracy as a function of marker selection density.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/*_Folds_*_Times_*.txt: Iteration-by-iteration predictive ability (Pearson r) and accuracy (r/h) logs for each model.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `data.table`, `snpReady`, `SNPRelate`, `rrBLUP`, `BGLR`, `dplyr`, `glmnet`, `randomForest`, `e1071`, `brnn`, `pedigree`, `reshape2`, `regress`, `MASS`

To install any missing packages in R:
```R
required_packages <- c("data.table", "snpReady", "SNPRelate", "rrBLUP", "BGLR", "dplyr", "glmnet", "randomForest", "e1071", "brnn", "pedigree", "reshape2", "regress", "MASS")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 03_genomic_prediction_cross_validation.R
```
or inside an R interactive session:
```R
source("03_genomic_prediction_cross_validation.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
