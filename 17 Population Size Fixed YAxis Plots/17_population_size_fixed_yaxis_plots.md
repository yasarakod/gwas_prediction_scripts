# Standardized Y-Axis Population Size Plotting

**Script File:** `17_population_size_fixed_yaxis_plots.R`  
**Location:** `17 Population Size Fixed YAxis Plots/17_population_size_fixed_yaxis_plots.R`

---

## 1. Brief Description
Generates population size comparison boxplots with standardized fixed y-axis scales ([-0.25, 0.65]) and bounded significance letter placement for uniform visual comparison across traits.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `15_population_size_friedman_analysis.R`
- **Module Folder:** `15 Population Size Friedman Analysis`
- **Role in Pipeline:** Generates Table_GBLUP_Times_10_*_top_1000_pop_folds_comparison.csv containing structured accuracy data across sample sizes and folds.

### 2. `09_gwas_markers_population_sizes_prediction.R (Alternative)`
- **Module Folder:** `09 GWAS Markers Population Sizes Prediction`
- **Role in Pipeline:** Provides raw sample size logs (GBLUP_Folds_*_top_1000.txt) if pre-compiled comparison tables are not available.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `Table_GBLUP_Times_10_*_top_1000_pop_folds_comparison.csv`
- **Data Description:** Consolidated population size comparison table.
- **Expected File Format:** Comma-Separated Values (CSV) with header.
- **Header & Structure:** Columns: Population_Size, Folds, Accuracy, Predictive_Ability.
- **Value Types & Encodings:** Numeric factors and metrics.
- **How the Code Processes It:** Read via read.csv(), letter positioning calculated with pmin(max_val + 0.02 * y_range, Y_MAX - 0.02), and rendered with coord_cartesian(ylim=c(-0.25, 0.65)).

### 2. `GBLUP_Folds_*_top_1000.txt`
- **Data Description:** Raw sample size iteration logs (optional fallback).
- **Expected File Format:** Single-column text files.
- **Header & Structure:** Iteration accuracy scores.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Fallback input if pre-compiled CSV table is not present.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/GBLUP_Boxplot_PredictiveAbility_<trait>_FixedY.png: Standardized y-axis predictive ability boxplot.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `dplyr`, `reshape2`, `PMCMRplus`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "dplyr", "reshape2", "PMCMRplus")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 17_population_size_fixed_yaxis_plots.R
```
or inside an R interactive session:
```R
source("17_population_size_fixed_yaxis_plots.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
