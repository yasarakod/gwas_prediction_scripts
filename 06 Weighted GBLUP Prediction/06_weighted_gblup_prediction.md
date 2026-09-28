# Iterative Weighted GBLUP (wGBLUP) Pipeline

**Script File:** `06_weighted_gblup_prediction.R`  
**Location:** `06 Weighted GBLUP Prediction/06_weighted_gblup_prediction.R`

---

## 1. Brief Description
Implements iterative Weighted GBLUP (wGBLUP) where SNP variance weights (d_i) are computed iteratively from marker effects to build trait-specific weighted G-matrices.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Produces genotype dosage matrix M, phenotype observations pheno, and marker map hapmap.

### 2. `GWAS Wald Association Pipeline (External / PLINK)`
- **Module Folder:** `data/`
- **Role in Pipeline:** Generates wald_DPC_date.txt, wald_DPC_time.txt, and wald_Is_Dead.txt containing single-SNP Wald test summary statistics used to calculate initial SNP weights.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Additive genotype dosage matrix (N individuals x M SNPs).
- **Expected File Format:** Numeric matrix object.
- **Header & Structure:** Rows = individual sample IDs, Columns = SNP marker IDs.
- **Value Types & Encodings:** 0, 1, 2 numeric dosages.
- **How the Code Processes It:** Used with diagonal weight matrix D to compute weighted genomic relationship matrix: G = (M %*% D %*% t(M)) / sum(2 * p * (1-p) * d).

### 2. `pheno`
- **Data Description:** Phenotypic records for target mortality traits.
- **Expected File Format:** Dataframe object.
- **Header & Structure:** Contains individual observations for DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Continuous and binary observations.
- **How the Code Processes It:** Response variable in mixed model equations.

### 3. `wald_DPC_date.txt / wald_DPC_time.txt / wald_Is_Dead.txt`
- **Data Description:** Genome-Wide Association Study (GWAS) Wald test summary statistics.
- **Expected File Format:** Tab-delimited text file with header.
- **Header & Structure:** Columns: SNP_ID, Chromosome, Position, Beta, SE, Wald_Statistic, P_Value.
- **Value Types & Encodings:** Floating-point test statistics and significance p-values.
- **How the Code Processes It:** Used to initialize SNP weights based on single-marker association significance.

### 4. `marker_weights_*.rds (e.g. marker_weights_DPC_time_top_1000.rds)`
- **Data Description:** Precalculated iterative SNP variance weights vector.
- **Expected File Format:** Binary RDS file containing named numeric vector.
- **Header & Structure:** Named vector of length equal to number of evaluated markers.
- **Value Types & Encodings:** Positive floating-point scaling weights.
- **How the Code Processes It:** Constructs the diagonal weight matrix D in iterative wGBLUP rounds.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/WGBLUP_Folds_*_Times_*.txt: Cross-validation prediction accuracy logs for weighted GBLUP.`
- `results/marker_weights_<trait>_<marker_size>.rds: Updated iterative SNP variance weights.`

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
Rscript 06_weighted_gblup_prediction.R
```
or inside an R interactive session:
```R
source("06_weighted_gblup_prediction.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
