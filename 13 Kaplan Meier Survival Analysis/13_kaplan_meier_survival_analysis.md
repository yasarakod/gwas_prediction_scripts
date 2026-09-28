# Kaplan-Meier Survival Analysis for Disease Challenge

**Script File:** `13_kaplan_meier_survival_analysis.R`  
**Location:** `13 Kaplan Meier Survival Analysis/13_kaplan_meier_survival_analysis.R`

---

## 1. Brief Description
Computes non-parametric Kaplan-Meier survival probability curves, median survival durations, and survival distributions for viral challenge mortality traits.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Cleans and prepares pheno (or can run directly on raw phenotypes.csv).

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pheno / pheno.rds / phenotypes.csv`
- **Data Description:** Phenotypic records containing time-to-event and censoring indicators.
- **Expected File Format:** Dataframe or CSV file with header.
- **Header & Structure:** Must contain columns: DPC_time (continuous time in hours), DPC_date (discrete time in days), and Is_Dead (event status indicator).
- **Value Types & Encodings:** Time: Positive continuous/integer values; Status: 1 = event (mortality), 0 = censored (survived to end of challenge).
- **How the Code Processes It:** Supplied to survival::Surv(time = DPC_time, event = Is_Dead) and fitted via survival::survfit().

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/kaplan_meier_dpc_time.png: Survival probability curve over hours post-challenge with confidence intervals.`
- `results/kaplan_meier_dpc_date.png: Survival probability curve over days post-challenge.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `survival`, `survminer`, `dplyr`

To install any missing packages in R:
```R
required_packages <- c("survival", "survminer", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 13_kaplan_meier_survival_analysis.R
```
or inside an R interactive session:
```R
source("13_kaplan_meier_survival_analysis.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
