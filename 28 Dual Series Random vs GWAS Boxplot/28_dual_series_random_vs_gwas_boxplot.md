# Dual-Series Random vs GWAS Marker Boxplot

**Script File:** `28_dual_series_random_vs_gwas_boxplot.R`  
**Location:** `28 Dual Series Random vs GWAS Boxplot/28_dual_series_random_vs_gwas_boxplot.R`

---

## 1. Brief Description
Creates dual-series color-coded boxplots directly contrasting Random marker sets (purple) vs GWAS-prioritized marker sets (green) across subset sizes.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `14_gwas_vs_random_marker_comparison.R`
- **Module Folder:** `14 GWAS vs Random Marker Comparison`
- **Role in Pipeline:** Processes raw marker density prediction runs to produce summary comparison tables (Table_GBLUP_Folds_5_Times_10_*_marker.csv / 1/2/3. Table_GBLUP_Folds_5_Times_10_*_marker.csv).

### 2. `03_genomic_prediction_cross_validation.R (Upstream)`
- **Module Folder:** `03 Genomic Prediction Cross Validation`
- **Role in Pipeline:** Generates the underlying cross-validation prediction iterations across Random and GWAS marker subsets.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `Table_GBLUP_Folds_5_Times_10_*_marker.csv (e.g. 1/2/3. Table_GBLUP_Folds_5_Times_10_*_marker.csv)`
- **Data Description:** Consolidated marker comparison data table containing prediction accuracy across marker subsets.
- **Expected File Format:** Comma-Separated Values (CSV) with header.
- **Header & Structure:** Rows = iterations / subset sizes, Columns = subset names (e.g. 50, top_50, 100, top_100, 500, top_500, 1000, top_1000, etc.).
- **Value Types & Encodings:** Numeric correlation values.
- **How the Code Processes It:** Reshaped to extract numeric subset sizes and category labels (Random vs GWAS) for ggplot2 dual-series boxplot rendering.

### 2. `Table_GBLUP_Folds_5_Times_10_*_marker_comparison.csv`
- **Data Description:** Alternative marker comparison summary table.
- **Expected File Format:** CSV table.
- **Header & Structure:** Marker density vs accuracy.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Alternative input source.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Ability_Plot_DPC_time_WGBLUP_vs_GBLUP_markers.png: Dual-series boxplot comparing Random vs GWAS markers.`
- `results/Dual_Series_Random_vs_GWAS_Markers.png: High-resolution dual-series comparison figure.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `dplyr`, `stringr`, `reshape2`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "dplyr", "stringr", "reshape2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 28_dual_series_random_vs_gwas_boxplot.R
```
or inside an R interactive session:
```R
source("28_dual_series_random_vs_gwas_boxplot.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
