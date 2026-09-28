# Enhanced Friedman & Nemenyi Ranking Analysis with CLD

**Script File:** `21_enhanced_friedman_ranking_analysis.R`  
**Location:** `21 Enhanced Friedman Ranking Analysis/21_enhanced_friedman_ranking_analysis.R`

---

## 1. Brief Description
Advanced Friedman-Nemenyi analysis generating Compact Letter Displays (CLD), Kendall W concordance coefficients, and detailed statistical text reports.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `07_vrma_prediction_models.R / 08_vrma_prediction_excel_export.R`
- **Module Folder:** `07 VRMA Prediction Models / 08 VRMA Prediction Excel Export`
- **Role in Pipeline:** Produces Table_Folds_5_Times_10_*.csv and raw iteration text files (*_Folds_5_Times_10_*.txt).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `Table_Folds_5_Times_10_*_top_40000.csv`
- **Data Description:** Consolidated cross-validation accuracy table across models.
- **Expected File Format:** CSV file with header.
- **Header & Structure:** Columns: Iteration, Method, Accuracy.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Analyzed to extract rank matrices, pairwise p-value matrices, and CLD letter groupings.

### 2. `*_Folds_5_Times_10_*.txt`
- **Data Description:** Raw method iteration text files.
- **Expected File Format:** Single-column text files.
- **Header & Structure:** Iteration accuracy values.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Fallback input files for matrix assembly.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Friedman_Nemenyi_Statistical_Report.txt: Detailed statistical report with Kendall W and pairwise comparisons.`
- `results/Friedman_<trait>_analysis.png: Enhanced significance boxplot with ranking letters.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `ggpubr`, `PMCMRplus`, `dplyr`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "ggpubr", "PMCMRplus", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 21_enhanced_friedman_ranking_analysis.R
```
or inside an R interactive session:
```R
source("21_enhanced_friedman_ranking_analysis.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
