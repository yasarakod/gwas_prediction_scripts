# GBLUP Population Size & Fold Trajectory Plots

**Script File:** `18_gblup_popsize_folds_trajectories.R`  
**Location:** `18 GBLUP PopSize Folds Trajectories/18_gblup_popsize_folds_trajectories.R`

---

## 1. Brief Description
Plots mean accuracy and predictive ability performance trajectories across increasing training sample sizes, comparing 3-fold vs 5-fold cross-validation schemes.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `09_gwas_markers_population_sizes_prediction.R`
- **Module Folder:** `09 GWAS Markers Population Sizes Prediction`
- **Role in Pipeline:** Generates cross-validation iteration logs across population sizes (GBLUP_Folds_3/5_Times_10_*).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `GBLUP_Folds_3/5_Times_10_*_top_1000.txt`
- **Data Description:** Cross-validation result logs across sample sizes for DPC_time, DPC_date, and Is_Dead.
- **Expected File Format:** Single-column numeric text files.
- **Header & Structure:** Iteration scores across 50 iterations per sample size.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Aggregated into mean and standard error metrics across sample sizes to construct trajectory line charts.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Table_GBLUP_Combined.csv: Combined accuracy dataframe across folds and sizes.`
- `results/Plot_GBLUP_Combined.png: Mean trajectory line plot for prediction accuracy.`
- `results/Ability_Plot_GBLUP_Combined.png: Mean trajectory line plot for predictive ability.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `dplyr`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 18_gblup_popsize_folds_trajectories.R
```
or inside an R interactive session:
```R
source("18_gblup_popsize_folds_trajectories.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
