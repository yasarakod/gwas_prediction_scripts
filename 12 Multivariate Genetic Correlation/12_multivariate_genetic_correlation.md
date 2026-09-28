# Multivariate Genetic & Phenotypic Correlation Analysis

**Script File:** `12_multivariate_genetic_correlation.R`  
**Location:** `12 Multivariate Genetic Correlation/12_multivariate_genetic_correlation.R`

---

## 1. Brief Description
Fits multivariate mixed models to estimate genetic, phenotypic, and environmental correlation matrices among mortality traits and morphological traits, generating graphical correlograms.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Produces M and multi-trait phenotype dataset pheno.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Produces A_matrix.rds for pedigree multivariate models.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pheno / pheno.rds`
- **Data Description:** Multi-trait phenotype dataframe.
- **Expected File Format:** Serialized dataframe or CSV.
- **Header & Structure:** N rows x P traits (DPC_time, DPC_date, Is_Dead, Length, Weight, Width).
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Supplied as multivariate response matrix Y to fit multi-trait mixed models.

### 2. `M`
- **Data Description:** Genotype dosage matrix.
- **Expected File Format:** Numeric matrix (N x M).
- **Header & Structure:** 0, 1, 2 allele dosages.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Used to construct genomic covariance matrix G for multi-trait GBLUP.

### 3. `A_matrix.rds`
- **Data Description:** Pedigree relationship matrix.
- **Expected File Format:** Binary RDS symmetric matrix.
- **Header & Structure:** N x N matrix.
- **Value Types & Encodings:** Kinship coefficients.
- **How the Code Processes It:** Covariance structure in multivariate pedigree models.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/genetic_correlations_all_methods.rds: Estimated genetic and phenotypic covariance and correlation matrices.`
- `results/genetic_correlations_gblup.png: Visual correlogram plot displaying genetic correlation coefficients.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `rrBLUP`, `BGLR`, `sommer`, `AGHmatrix`, `corrplot`, `ggplot2`

To install any missing packages in R:
```R
required_packages <- c("rrBLUP", "BGLR", "sommer", "AGHmatrix", "corrplot", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 12_multivariate_genetic_correlation.R
```
or inside an R interactive session:
```R
source("12_multivariate_genetic_correlation.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
