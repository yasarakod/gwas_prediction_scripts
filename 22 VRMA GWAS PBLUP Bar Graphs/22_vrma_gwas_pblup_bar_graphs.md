# Grouped Bar Chart Visualization with Friedman Statistics

**Script File:** `22_vrma_gwas_pblup_bar_graphs.R`  
**Location:** `22 VRMA GWAS PBLUP Bar Graphs/22_vrma_gwas_pblup_bar_graphs.R`

---

## 1. Brief Description
Generates grouped bar charts comparing prediction methods (including PBLUP) across viral challenge traits, displaying error bars (mean +- SD) and Kendall W concordance in legend.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `07_vrma_prediction_models.R / 08_vrma_prediction_excel_export.R`
- **Module Folder:** `07 VRMA Prediction Models / 08 VRMA Prediction Excel Export`
- **Role in Pipeline:** Executes 10-model cross-validation with top 1000 GWAS markers and generates iteration logs (*_Folds_5_Times_10_*_top_1000.txt) and Prediction_Accuracy_Iterations_AllTraits.csv.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `*_Folds_5_Times_10_*_top_1000.txt`
- **Data Description:** Iteration accuracy logs across 10 prediction methods (PBLUP, GBLUP, BayesA/B/C, BRR, BL, EN, RF, RR) for 1000 GWAS markers.
- **Expected File Format:** Single-column text files.
- **Header & Structure:** 50 numeric accuracy scores per file.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Loaded to compute method means, standard deviations, and Kendall W concordance coefficients.

### 2. `Prediction_Accuracy_Iterations_AllTraits.csv / Method_Performance_Summary.csv`
- **Data Description:** Consolidated iteration-level accuracy table across all methods and traits.
- **Expected File Format:** CSV file with header.
- **Header & Structure:** Columns: Method, Trait, Accuracy.
- **Value Types & Encodings:** Numeric metrics.
- **How the Code Processes It:** Alternative consolidated source for summary statistics.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Grouped_Bar_Chart_Main.png: Main grouped bar chart comparing methods across traits.`
- `results/Grouped_Bar_Chart_With_Stats.png: Statistical grouped bar chart with Kendall concordance in subtitle.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `PMCMRplus`, `dplyr`, `rcompanion`, `gridExtra`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "PMCMRplus", "dplyr", "rcompanion", "gridExtra")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 22_vrma_gwas_pblup_bar_graphs.R
```
or inside an R interactive session:
```R
source("22_vrma_gwas_pblup_bar_graphs.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
