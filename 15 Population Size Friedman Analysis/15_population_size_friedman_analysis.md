# Training Population Size Performance & Friedman Testing

**Script File:** `15_population_size_friedman_analysis.R`  
**Location:** `15 Population Size Friedman Analysis/15_population_size_friedman_analysis.R`

---

## 1. Brief Description
Evaluates prediction performance across training sample sizes (50 to 474 individuals) and cross-validation schemes (3-fold vs 5-fold) using Friedman ANOVA and Nemenyi ranking letter tests.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `09_gwas_markers_population_sizes_prediction.R`
- **Module Folder:** `09 GWAS Markers Population Sizes Prediction`
- **Role in Pipeline:** Runs cross-validation across sample sizes (50 to 474) with 1000 GWAS markers to generate GBLUP_Folds_3/5_Times_10_*_<popsize>_top_1000.txt logs.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `GBLUP_Folds_3_Times_10_*_top_1000.txt & GBLUP_Folds_5_Times_10_*_top_1000.txt`
- **Data Description:** Cross-validation iteration logs across training population sizes (50, 100, 150, 200, 250, 300, 350, 400, 474 individuals).
- **Expected File Format:** Single-column text files.
- **Header & Structure:** Rows contain numeric accuracy scores per cross-validation iteration.
- **Value Types & Encodings:** Floating-point values.
- **How the Code Processes It:** Constructs population size comparison matrices for statistical ranking.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Table_<pmethod><namef>.csv: Summary table of sample size performance.`
- `results/Friedman_Results_<pmethod><namef>.csv: Friedman ANOVA statistics.`
- `results/PopFolds_Rankings_<pmethod><namef>.csv: Significance ranking letters across sample sizes.`
- `results/Accuracy_Plot2_<pmethod><namef>.png: Accuracy multi-panel chart.`
- `results/Ability_Plot2_<pmethod><namef>.png: Predictive ability chart.`

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
Rscript 15_population_size_friedman_analysis.R
```
or inside an R interactive session:
```R
source("15_population_size_friedman_analysis.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
