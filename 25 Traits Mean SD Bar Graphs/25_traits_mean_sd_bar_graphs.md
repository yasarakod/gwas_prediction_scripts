# Mean and SD Accuracy Bar Charts for Top 1000 Mb Markers

**Script File:** `25_traits_mean_sd_bar_graphs.R`  
**Location:** `25 Traits Mean SD Bar Graphs/25_traits_mean_sd_bar_graphs.R`

---

## 1. Brief Description
Calculates mean and standard deviation of predictive ability and accuracy across challenge traits with top 1000 Mb GWAS markers and produces grouped bar graphs.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates M, pheno, and top_1000_snps.txt.

### 2. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Generates A_matrix.rds.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Genotype dosage matrix.
- **Expected File Format:** Numeric matrix object.
- **Header & Structure:** N individuals x M SNPs.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Model feature matrix.

### 2. `pheno`
- **Data Description:** Phenotypic records.
- **Expected File Format:** Dataframe.
- **Header & Structure:** Columns: DPC_time, DPC_date, Is_Dead.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Model response variables.

### 3. `A_matrix.rds`
- **Data Description:** Additive relationship matrix from pedigree.
- **Expected File Format:** Binary RDS matrix.
- **Header & Structure:** N x N matrix.
- **Value Types & Encodings:** Kinship coefficients.
- **How the Code Processes It:** PBLUP covariance.

### 4. `top_1000_snps.txt`
- **Data Description:** Top 1000 marker list.
- **Expected File Format:** Single-column text file.
- **Header & Structure:** 1,000 SNP names.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Marker filtering.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Bargraph_Traits_Mean_SD_Top1000Mb.png: Mean and SD summary bar chart.`

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
Rscript 25_traits_mean_sd_bar_graphs.R
```
or inside an R interactive session:
```R
source("25_traits_mean_sd_bar_graphs.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
