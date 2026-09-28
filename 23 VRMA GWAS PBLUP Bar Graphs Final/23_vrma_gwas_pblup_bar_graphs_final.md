# Polished Publication-Ready Grouped Bar Charts

**Script File:** `23_vrma_gwas_pblup_bar_graphs_final.R`  
**Location:** `23 VRMA GWAS PBLUP Bar Graphs Final/23_vrma_gwas_pblup_bar_graphs_final.R`

---

## 1. Brief Description
Generates high-resolution publication-quality grouped bar charts with Viridis color palettes, customized theme aesthetics, and error bars.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `07_vrma_prediction_models.R / 08_vrma_prediction_excel_export.R`
- **Module Folder:** `07 VRMA Prediction Models / 08 VRMA Prediction Excel Export`
- **Role in Pipeline:** Generates the multi-model cross-validation iteration logs (*_Folds_5_Times_10_*_top_1000.txt) and Prediction_Accuracy_Iterations_AllTraits.csv.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `*_Folds_5_Times_10_*_top_1000.txt`
- **Data Description:** Iteration accuracy logs for 1000 GWAS markers across 10 prediction methods.
- **Expected File Format:** Single-column numeric text files.
- **Header & Structure:** 50 values per file.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Aggregated into summary data frame with mean and standard deviation.

### 2. `Prediction_Accuracy_Iterations_AllTraits.csv`
- **Data Description:** Consolidated prediction iteration table.
- **Expected File Format:** CSV table with header.
- **Header & Structure:** Columns: Method, Trait, Accuracy.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Pre-compiled input table.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Final_Grouped_Bar_Chart_Legend_Spacing_Viridis.png: 300 DPI publication bar chart in Viridis color palette.`
- `results/Publication_Ready_Better_Spacing.png: High-resolution figure with adjusted label spacing.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `PMCMRplus`, `dplyr`, `viridis`, `rcompanion`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "PMCMRplus", "dplyr", "viridis", "rcompanion")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 23_vrma_gwas_pblup_bar_graphs_final.R
```
or inside an R interactive session:
```R
source("23_vrma_gwas_pblup_bar_graphs_final.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
