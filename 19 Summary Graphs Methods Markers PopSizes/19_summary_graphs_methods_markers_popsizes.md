# Multi-Factor Summary Graphs & Comparison Tables

**Script File:** `19_summary_graphs_methods_markers_popsizes.R`  
**Location:** `19 Summary Graphs Methods Markers PopSizes/19_summary_graphs_methods_markers_popsizes.R`

---

## 1. Brief Description
Synthesizes raw cross-validation results across prediction methods, marker densities, and sample sizes into multi-panel graphical summaries and comparative tables.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `03_genomic_prediction_cross_validation.R / 07_vrma_prediction_models.R`
- **Module Folder:** `03 Genomic Prediction Cross Validation / 07 VRMA Prediction Models`
- **Role in Pipeline:** Generates raw cross-validation log files (<Method>_Folds_*_Times_*.txt) across prediction models and marker subsets.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `<Method>_Folds_*_Times_*.txt`
- **Data Description:** Cross-validation output files across all evaluated algorithms (GBLUP, PBLUP, BayesA/B/C, BL, BRR, EN, RF, RR, WGBLUP).
- **Expected File Format:** Single-column text files.
- **Header & Structure:** Iteration accuracy values.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Read via read.table() and compiled into summary tables and multi-panel figures.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Table<namef>.csv: Tabular summary of model comparisons.`
- `results/Plot<namef>.png: Accuracy multi-factor comparison plot.`
- `results/Ability_Plot<namef>.png: Predictive ability comparison plot.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `argparse`, `ggplot2`, `reshape2`, `ggpubr`

To install any missing packages in R:
```R
required_packages <- c("argparse", "ggplot2", "reshape2", "ggpubr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 19_summary_graphs_methods_markers_popsizes.R
```
or inside an R interactive session:
```R
source("19_summary_graphs_methods_markers_popsizes.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
