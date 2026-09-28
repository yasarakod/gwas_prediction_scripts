# Command-Line Modular Genomic Prediction Pipeline

**Script File:** `05_genomic_prediction_cli_module.R`  
**Location:** `05 Genomic Prediction CLI Module/05_genomic_prediction_cli_module.R`

---

## 1. Brief Description
Modular command-line script leveraging argparse to execute parametric genomic prediction cross-validation runs with customizable flags for traits, folds, iterations, and markers.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates the genotype dosage matrix M and phenotype data file pheno.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Generates completed_pedigree.rds and A_matrix.rds.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pheno / pheno.rds`
- **Data Description:** Target phenotype dataframe for model fitting.
- **Expected File Format:** RDS serialized dataframe or tabular text file.
- **Header & Structure:** Rows = individuals, Columns = evaluated traits.
- **Value Types & Encodings:** Numeric phenotypic observations.
- **How the Code Processes It:** Specified via --pheno_file CLI argument.

### 2. `M`
- **Data Description:** Genotype dosage matrix (N individuals x M SNPs).
- **Expected File Format:** RDS serialized matrix or text file.
- **Header & Structure:** Dosage matrix with individual IDs as row names.
- **Value Types & Encodings:** Numeric 0, 1, 2.
- **How the Code Processes It:** Specified via --geno_file CLI argument.

### 3. `A_matrix.rds`
- **Data Description:** Numerator additive pedigree relationship matrix.
- **Expected File Format:** Binary RDS file.
- **Header & Structure:** N x N square covariance matrix.
- **Value Types & Encodings:** Numeric kinship coefficients.
- **How the Code Processes It:** Loaded for pedigree-based BLUP models when requested.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/<Method>_Folds_<fold>_Times_<iter>_<trait>_<marker_size>.txt: Output logs containing cross-validation accuracy values.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `argparse`, `rrBLUP`, `BGLR`, `glmnet`, `randomForest`, `e1071`, `brnn`, `AGHmatrix`

To install any missing packages in R:
```R
required_packages <- c("argparse", "rrBLUP", "BGLR", "glmnet", "randomForest", "e1071", "brnn", "AGHmatrix")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 05_genomic_prediction_cli_module.R
```
or inside an R interactive session:
```R
source("05_genomic_prediction_cli_module.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
