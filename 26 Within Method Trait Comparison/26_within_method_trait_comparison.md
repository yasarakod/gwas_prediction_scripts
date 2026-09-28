# Within-Method Trait Predictability Comparison

**Script File:** `26_within_method_trait_comparison.R`  
**Location:** `26 Within Method Trait Comparison/26_within_method_trait_comparison.R`

---

## 1. Brief Description
Statistically evaluates whether certain traits are more predictable than others within each prediction algorithm, generating performance heatmaps and trait ranking letter tables.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `07_vrma_prediction_models.R / 08_vrma_prediction_excel_export.R`
- **Module Folder:** `07 VRMA Prediction Models / 08 VRMA Prediction Excel Export`
- **Role in Pipeline:** Runs cross-validation across all 10 prediction algorithms for DPC_time, DPC_date, and Is_Dead to generate the required iteration logs (*_Folds_5_Times_10_*_top_1000.txt) and Prediction_Accuracy_Iterations_AllTraits.csv.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `*_Folds_5_Times_10_*_top_1000.txt`
- **Data Description:** Iteration accuracy logs across traits (DPC_time, DPC_date, Is_Dead) for each model.
- **Expected File Format:** Single-column text files.
- **Header & Structure:** 50 accuracy values per model-trait pair.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Evaluated within each model via Friedman test and manual/multcompView Nemenyi post-hoc grouping.

### 2. `Prediction_Accuracy_Iterations_AllTraits.csv`
- **Data Description:** Consolidated iteration table.
- **Expected File Format:** CSV table with header.
- **Header & Structure:** Columns: Method, Trait, Accuracy.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Consolidated input source.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Predictive_Ability_With_Significance_Ranked.png: Ranked comparison chart across traits.`
- `results/Trait_Performance_Heatmap_Ranked.png: Trait vs Method performance heatmap.`
- `results/Method_Performance_Summary.csv: Comprehensive statistical summary table with ranking letters.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `PMCMRplus`, `dplyr`, `tidyr`, `multcompView`, `viridis`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "PMCMRplus", "dplyr", "tidyr", "multcompView", "viridis")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 26_within_method_trait_comparison.R
```
or inside an R interactive session:
```R
source("26_within_method_trait_comparison.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
