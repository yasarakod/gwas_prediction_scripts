# Population Size Impact on GWAS Marker Genomic Prediction

**Script File:** `09_gwas_markers_population_sizes_prediction.R`  
**Location:** `09 GWAS Markers Population Sizes Prediction/09_gwas_markers_population_sizes_prediction.R`

---

## 1. Brief Description
Assesses genomic prediction accuracy as a function of training population sample size (50, 100, 150, 200, 250, 300, 350, 400, 474 individuals) using top 1000 GWAS markers.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates M, pheno, and top_1000_snps.txt.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Genotype dosage matrix.
- **Expected File Format:** Numeric matrix object.
- **Header & Structure:** Rows = individuals, Columns = SNPs.
- **Value Types & Encodings:** 0, 1, 2 allele counts.
- **How the Code Processes It:** Subsampled to specified training population sizes in each cross-validation iteration.

### 2. `pheno`
- **Data Description:** Phenotypic records.
- **Expected File Format:** Dataframe.
- **Header & Structure:** Contains response traits (DPC_time, DPC_date, Is_Dead).
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Subsampled alongside genotype matrix M.

### 3. `top_1000_snps.txt`
- **Data Description:** List of top 1000 GWAS markers.
- **Expected File Format:** Plain text list.
- **Header & Structure:** 1,000 lines of SNP IDs.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Restricts genomic prediction to top GWAS markers.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/GBLUP_Folds_*_Times_10_*_<popsize>_top_1000.txt: Predictive accuracy logs for each evaluated training sample size.`

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
Rscript 09_gwas_markers_population_sizes_prediction.R
```
or inside an R interactive session:
```R
source("09_gwas_markers_population_sizes_prediction.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
