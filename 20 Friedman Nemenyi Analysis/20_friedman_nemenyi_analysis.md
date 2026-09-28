# Friedman & Nemenyi Non-Parametric Model Ranking

**Script File:** `20_friedman_nemenyi_analysis.R`  
**Location:** `20 Friedman Nemenyi Analysis/20_friedman_nemenyi_analysis.R`

---

## 1. Brief Description
Performs non-parametric Friedman two-way ANOVA by ranks and Nemenyi post-hoc pairwise comparisons to rank prediction methods and establish statistical significance groups.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `07_vrma_prediction_models.R / 08_vrma_prediction_excel_export.R`
- **Module Folder:** `07 VRMA Prediction Models / 08 VRMA Prediction Excel Export`
- **Role in Pipeline:** Generates iteration-level prediction accuracy results and summary tables (Table_Folds_5_Times_10_*_top_40000.csv / Table_Folds_5_Times_10_*_top_1000.csv).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `Table_Folds_5_Times_10_*_top_40000.csv / Table_Folds_5_Times_10_*_top_1000.csv`
- **Data Description:** Summary tables containing iteration-level prediction accuracy scores across methods.
- **Expected File Format:** Comma-Separated Values (CSV) with header.
- **Header & Structure:** Columns: Iteration, Method, Accuracy, Trait.
- **Value Types & Encodings:** Numeric and factor columns.
- **How the Code Processes It:** Reshaped into a block matrix (methods as columns, iterations as blocks) and passed to PMCMRplus::frdAllPairsNemenyiTest().

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Friedman_Results_*.csv: Friedman test chi-squared statistic and p-value.`
- `results/Method_Rankings_*.csv: Mean ranks and significance grouping letters.`
- `results/Plot_with_Rankings_*.png: Boxplot annotated with method rankings.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `tidyverse`, `PMCMRplus`, `rstatix`, `ggplot2`

To install any missing packages in R:
```R
required_packages <- c("tidyverse", "PMCMRplus", "rstatix", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 20_friedman_nemenyi_analysis.R
```
or inside an R interactive session:
```R
source("20_friedman_nemenyi_analysis.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
