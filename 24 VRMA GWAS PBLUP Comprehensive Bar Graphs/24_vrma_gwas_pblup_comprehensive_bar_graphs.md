# End-to-End VRMA Prediction & Bar Graph Generation

**Script File:** `24_vrma_gwas_pblup_comprehensive_bar_graphs.R`  
**Location:** `24 VRMA GWAS PBLUP Comprehensive Bar Graphs/24_vrma_gwas_pblup_comprehensive_bar_graphs.R`

---

## 1. Brief Description
Integrated workflow combining data preprocessing, 10-model cross-validation execution for 1000 GWAS markers and PBLUP, and automated summary bar chart creation.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Produces genotype dosage matrix M, phenotype dataframe pheno, and top_1000_snps.txt.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Produces A_matrix.rds and completed_pedigree.rds.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Genotype dosage matrix (N individuals x M SNPs).
- **Expected File Format:** Numeric matrix object.
- **Header & Structure:** 0, 1, 2 allele dosages.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Feature matrix in cross-validation.

### 2. `pheno`
- **Data Description:** Cleaned phenotype observations.
- **Expected File Format:** Dataframe.
- **Header & Structure:** Columns: DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Response variables in model training.

### 3. `A_matrix.rds / completed_pedigree.rds`
- **Data Description:** Additive relationship matrix and sorted pedigree dataset.
- **Expected File Format:** Binary RDS objects.
- **Header & Structure:** Covariance matrix and parentage table.
- **Value Types & Encodings:** Numeric / character.
- **How the Code Processes It:** PBLUP covariance structure.

### 4. `top_1000_snps.txt`
- **Data Description:** Top 1000 GWAS marker list.
- **Expected File Format:** Single-column text file.
- **Header & Structure:** 1,000 SNP names.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Filters genotype matrix to top 1000 markers.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/*_Folds_5_Times_10_*_top_1000.txt: Raw cross-validation output logs.`
- `results/Bargraph_Traits_Mean_SD_Top1000_PBLUP.png: Summary bar chart of mean prediction accuracy.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `data.table`, `snpReady`, `SNPRelate`, `rrBLUP`, `BGLR`, `dplyr`, `glmnet`, `randomForest`, `e1071`, `brnn`, `pedigree`, `reshape2`, `regress`, `MASS`, `ggplot2`

To install any missing packages in R:
```R
required_packages <- c("data.table", "snpReady", "SNPRelate", "rrBLUP", "BGLR", "dplyr", "glmnet", "randomForest", "e1071", "brnn", "pedigree", "reshape2", "regress", "MASS", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 24_vrma_gwas_pblup_comprehensive_bar_graphs.R
```
or inside an R interactive session:
```R
source("24_vrma_gwas_pblup_comprehensive_bar_graphs.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
