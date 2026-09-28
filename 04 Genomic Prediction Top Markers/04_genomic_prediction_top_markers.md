# Genomic Prediction with Top GWAS-Selected Markers

**Script File:** `04_genomic_prediction_top_markers.R`  
**Location:** `04 Genomic Prediction Top Markers/04_genomic_prediction_top_markers.R`

---

## 1. Brief Description
Performs 10-fold cross-validation across all prediction models utilizing top GWAS-prioritized SNP subsets (e.g. top 1000 markers vs ALL markers) for mortality challenge traits.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Required to produce the genotype dosage matrix M, clean phenotype dataset pheno, and top SNP lists.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Required to compute A_matrix.rds for pedigree-based model comparisons.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Numerical genotype dosage matrix (individuals x markers).
- **Expected File Format:** Preprocessed numeric matrix object.
- **Header & Structure:** N individuals x M SNPs. Values represent minor allele dosage (0, 1, 2).
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Filtered by top marker lists to construct training/test marker feature matrices.

### 2. `pheno / pheno.rds`
- **Data Description:** Phenotypic data frame for viral mortality traits.
- **Expected File Format:** Serialized dataframe.
- **Header & Structure:** Columns: DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Continuous and binary numeric variables.
- **How the Code Processes It:** Target response in 10-fold cross-validation.

### 3. `top_1000_snps.txt`
- **Data Description:** List of the top 1,000 most significant GWAS SNP markers from Wald test rankings.
- **Expected File Format:** Single-column text file.
- **Header & Structure:** 1,000 lines, each containing a single SNP name.
- **Value Types & Encodings:** Character strings matching matrix column names.
- **How the Code Processes It:** Filters genotype matrix columns: M_sub <- M[, top_1000_snps].

### 4. `top_ALL_snps.txt`
- **Data Description:** Complete inventory of all quality-controlled SNP markers across the genome.
- **Expected File Format:** Single-column text file.
- **Header & Structure:** Lines containing all QC-filtered marker names.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Benchmark baseline representing whole-genome marker prediction.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/*_Folds_5_Times_10_*_top_1000.txt: Cross-validation prediction logs for top 1000 GWAS markers across all 10 algorithms.`

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
Rscript 04_genomic_prediction_top_markers.R
```
or inside an R interactive session:
```R
source("04_genomic_prediction_top_markers.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
