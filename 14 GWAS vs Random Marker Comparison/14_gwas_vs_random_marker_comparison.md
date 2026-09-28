# GWAS vs Random Marker Statistical Comparison & Plotting

**Script File:** `14_gwas_vs_random_marker_comparison.R`  
**Location:** `14 GWAS vs Random Marker Comparison/14_gwas_vs_random_marker_comparison.R`

---

## 1. Brief Description
Performs statistical evaluations (Friedman non-parametric test and Nemenyi post-hoc ranking letter assignments) comparing prediction accuracy between GWAS-prioritized and Random marker density subsets.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `03_genomic_prediction_cross_validation.R`
- **Module Folder:** `03 Genomic Prediction Cross Validation`
- **Role in Pipeline:** Executes cross-validation across the marker density subsets (50 to ALL) for both Random and GWAS-prioritized markers to produce the raw GBLUP_Folds_5_Times_10_* log files.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `GBLUP_Folds_5_Times_10_*_top_*.txt & GBLUP_Folds_5_Times_10_*_*.txt`
- **Data Description:** Cross-validation accuracy logs generated across marker subsets (50, 100, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000, ALL) for GWAS vs Random markers.
- **Expected File Format:** Single-column plain text files (.txt).
- **Header & Structure:** 50 rows of numeric prediction correlations (5 folds x 10 repeats).
- **Value Types & Encodings:** Floating-point correlation values in [-1, 1].
- **How the Code Processes It:** Assembled into comparison matrix (rows = iterations, columns = marker subsets) to execute Friedman tests and Nemenyi ranking letter assignments.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Table_<pmethod><namef>.csv: Statistical summary table of mean accuracy and ranking.`
- `results/Friedman_Results_<pmethod><namef>.csv: Friedman test statistics, degrees of freedom, and p-values.`
- `results/Marker_Rankings_<pmethod><namef>.csv: Significance grouping letters (e.g. a, b, c).`
- `results/Accuracy_Plot_with_Rankings_*.png: Boxplot with annotated ranking letters.`
- `results/Ability_Plot_with_Rankings_*.png: Predictive ability boxplot with significance groupings.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `PMCMRplus`, `dplyr`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "PMCMRplus", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 14_gwas_vs_random_marker_comparison.R
```
or inside an R interactive session:
```R
source("14_gwas_vs_random_marker_comparison.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
