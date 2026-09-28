# Final Integrated VRMA Genomic & Pedigree Prediction Pipeline

**Script File:** `07_vrma_prediction_models.R`  
**Location:** `07 VRMA Prediction Models/07_vrma_prediction_models.R`

---

## 1. Brief Description
Comprehensive benchmarking pipeline integrating pedigree (PBLUP), genomic (GBLUP), Bayesian regression (BayesA/B/C, BL, BRR), and machine learning (RF, EN) on viral challenge datasets.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates the genotype dosage matrix M, phenotype dataset pheno, and top SNP lists.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Generates completed_pedigree.rds and additive relationship matrix A_matrix.rds.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Numeric genotype dosage matrix.
- **Expected File Format:** Matrix object (N x M).
- **Header & Structure:** Dosages of 0, 1, 2 for all individuals across markers.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Core feature matrix for all genomic prediction algorithms.

### 2. `pheno`
- **Data Description:** Phenotypic records for viral refractory challenge.
- **Expected File Format:** Dataframe with aligned sample rows.
- **Header & Structure:** Columns include DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Continuous and binary values.
- **How the Code Processes It:** Target response in 5-fold cross-validation with 10 replicates.

### 3. `A_matrix.rds`
- **Data Description:** Additive numerator relationship matrix derived from pedigree.
- **Expected File Format:** Binary RDS symmetric matrix.
- **Header & Structure:** N x N matrix with sample IDs.
- **Value Types & Encodings:** Floating-point additive kinship.
- **How the Code Processes It:** Supplied to PBLUP mixed model solver.

### 4. `top_1000_snps.txt / top_40000_snps.txt`
- **Data Description:** Ranked marker subset identifier lists.
- **Expected File Format:** Plain text files with one SNP ID per line.
- **Header & Structure:** Single column without header.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Subsets genotype matrix M before running model evaluations.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/*_Folds_5_Times_10_*.txt: Validation accuracy logs for all evaluated models and traits.`

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
Rscript 07_vrma_prediction_models.R
```
or inside an R interactive session:
```R
source("07_vrma_prediction_models.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
