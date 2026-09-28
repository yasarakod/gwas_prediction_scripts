# Advanced Population Size Statistical Comparison with Significance Letters

**Script File:** `16_population_size_statistical_analysis.R`  
**Location:** `16 Population Size Statistical Analysis/16_population_size_statistical_analysis.R`

---

## 1. Brief Description
Advanced sample size benchmarking with strict statistical significance thresholds (alpha = 0.01) for Nemenyi post-hoc grouping and manuscript figure generation.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `09_gwas_markers_population_sizes_prediction.R`
- **Module Folder:** `09 GWAS Markers Population Sizes Prediction`
- **Role in Pipeline:** Generates sample size cross-validation logs (GBLUP_Folds_3/5_Times_10_*_top_1000.txt).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `GBLUP_Folds_3/5_Times_10_*_top_1000.txt`
- **Data Description:** Accuracy logs across sample sizes for traits DPC_time, DPC_date, and Is_Dead.
- **Expected File Format:** Single-column numeric text files.
- **Header & Structure:** Iteration accuracy values across 3-fold and 5-fold CV.
- **Value Types & Encodings:** Floating-point correlation values.
- **How the Code Processes It:** Evaluated with alpha = 0.01 threshold to generate conservative letter groupings.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Table_GBLUP_Times_10_*_for_manuscript.csv: Manuscript summary statistics table.`
- `results/Friedman_Results_GBLUP_*_for_manuscript.csv: Statistical test results.`
- `results/PopSize_Rankings_GBLUP_*_for_manuscript.csv: Letter assignments at alpha = 0.01.`
- `results/Accuracy_PopSize_Comparison_*.png: Manuscript accuracy figure.`
- `results/Ability_PopSize_Comparison_*.png: Manuscript predictive ability figure.`

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
Rscript 16_population_size_statistical_analysis.R
```
or inside an R interactive session:
```R
source("16_population_size_statistical_analysis.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
