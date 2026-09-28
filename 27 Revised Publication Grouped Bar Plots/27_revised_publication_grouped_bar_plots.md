# Publication Grouped Bar Plots with Kendall Concordance

**Script File:** `27_revised_publication_grouped_bar_plots.R`  
**Location:** `27 Revised Publication Grouped Bar Plots/27_revised_publication_grouped_bar_plots.R`

---

## 1. Brief Description
Produces high-resolution grouped bar charts with significance grouping letters, error bars (mean +- SD), Kendall W concordance metrics, and dual PNG/PDF exports.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `08_vrma_prediction_excel_export.R / 26_within_method_trait_comparison.R`
- **Module Folder:** `08 VRMA Prediction Excel Export / 26 Within Method Trait Comparison`
- **Role in Pipeline:** Generates summary performance tables (Method_Performance_Summary.csv, Prediction_Accuracy_Iterations_AllTraits.csv, Table_Folds_5_Times_10_*.csv).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `Table_GBLUP_Combined.csv / Table_Folds_5_Times_10_*.csv`
- **Data Description:** Cross-validation accuracy summary tables.
- **Expected File Format:** CSV files with header.
- **Header & Structure:** Columns: Method, Trait, Accuracy.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Parsed to calculate method means, standard deviations, and significance letters.

### 2. `Method_Performance_Summary.csv / Prediction_Accuracy_Iterations_AllTraits.csv`
- **Data Description:** Performance summary and iteration table.
- **Expected File Format:** CSV files.
- **Header & Structure:** Contains summarized metrics and letter rankings.
- **Value Types & Encodings:** Numeric and character.
- **How the Code Processes It:** Supplied to ggplot2 rendering pipeline.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/Publication_Ready_Grouped_Bar_Chart.png: 300 DPI publication-ready figure.`
- `results/Publication_Ready_With_Stats.png: Grouped bar chart annotated with significance letters and Kendall W.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `ggplot2`, `reshape2`, `PMCMRplus`, `dplyr`, `viridis`

To install any missing packages in R:
```R
required_packages <- c("ggplot2", "reshape2", "PMCMRplus", "dplyr", "viridis")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 27_revised_publication_grouped_bar_plots.R
```
or inside an R interactive session:
```R
source("27_revised_publication_grouped_bar_plots.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
