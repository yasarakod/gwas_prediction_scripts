# Pedigree Additive Relationship Matrix Calculation

**Script File:** `01_pedigree_matrix_calculation.R`  
**Location:** `01 Pedigree Matrix Calculation/01_pedigree_matrix_calculation.R`

---

## 1. Brief Description
Calculates the numerator additive relationship matrix (A-matrix) and its mathematical inverse (Ainv) from historical pedigree records using the AGHmatrix package.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `None (Foundational Script)`
- **Role in Pipeline:** This is a foundational data generation script that operates directly on raw pedigree text/CSV files (e.g. pedigree.txt, Pedigree.csv). No prior R scripts need to be run.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pedigree.txt`
- **Data Description:** Three-column pedigree records detailing individual animal/plant identification and parentage (sire and dam).
- **Expected File Format:** Tab-delimited text file (.txt) or space-delimited table.
- **Header & Structure:** Columns: ID, Sire, Dam. First row contains header column names.
- **Value Types & Encodings:** Character IDs. Unknown sires or dams must be encoded as 0 or NA.
- **How the Code Processes It:** Read via read.table(..., header=TRUE, stringsAsFactors=FALSE) and supplied to AGHmatrix::Amatrix(pedigree, ploidy=2) to compute the dense N x N symmetric additive covariance matrix.

### 2. `Pedigree.csv`
- **Data Description:** Alternative comma-separated pedigree dataset containing identical parentage columns.
- **Expected File Format:** Comma-Separated Values (CSV) with header.
- **Header & Structure:** Columns: ID, Sire, Dam. Each row represents a distinct individual.
- **Value Types & Encodings:** Alphanumeric sample identifiers. Missing parental lineage encoded as NA or blank.
- **How the Code Processes It:** Serves as an interchangeable input format for pedigree data loading.

### 3. `Pedigree.xlsx`
- **Data Description:** Microsoft Excel workbook containing pedigree lineage records.
- **Expected File Format:** Excel spreadsheet (.xlsx) workbook.
- **Header & Structure:** Sheet 1 contains standard columns: ID, Sire, Dam.
- **Value Types & Encodings:** Standard spreadsheet character cells.
- **How the Code Processes It:** Can be loaded using readxl::read_excel() for automated data ingestion.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/completed_pedigree.rds: Formatted pedigree data frame with sorted, verified ancestral IDs.`
- `results/A_matrix.rds: N x N dense numeric additive relationship covariance matrix.`
- `results/Ainv_matrix.rds: Inverted relationship matrix (A^-1) used in mixed model equations.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `AGHmatrix`, `rrBLUP`

To install any missing packages in R:
```R
required_packages <- c("AGHmatrix", "rrBLUP")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 01_pedigree_matrix_calculation.R
```
or inside an R interactive session:
```R
source("01_pedigree_matrix_calculation.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
