# PBLUP vs GBLUP Heritability Comparison

**Script File:** `11_pblup_vs_gblup_heritability.R`  
**Location:** `11 PBLUP vs GBLUP Heritability/11_pblup_vs_gblup_heritability.R`

---

## 1. Brief Description
Calculates pedigree-based narrow-sense heritability using the numerator relationship matrix (A) and performs comparative evaluations against genomic heritability (GBLUP).

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `01_pedigree_matrix_calculation.R`
- **Module Folder:** `01 Pedigree Matrix Calculation`
- **Role in Pipeline:** Computes A_matrix.rds and completed_pedigree.rds.

### 2. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates M and pheno.

### 3. `10_gblup_heritability_estimation.R`
- **Module Folder:** `10 GBLUP Heritability Estimation`
- **Role in Pipeline:** Produces GBLUP heritability estimates (heritability_<trait>.rds) for direct comparison with PBLUP.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pheno / pheno.rds`
- **Data Description:** Phenotypic records for target traits.
- **Expected File Format:** Dataframe or serialized RDS.
- **Header & Structure:** Rows match relationship matrix dimensions; columns include DPC_time, DPC_date, Is_Dead, Length, Weight, Width.
- **Value Types & Encodings:** Numeric.
- **How the Code Processes It:** Response variable y in mixed.solve(y = y, K = A) and mixed.solve(y = y, K = G).

### 2. `A_matrix.rds`
- **Data Description:** Numerator additive pedigree relationship matrix (A).
- **Expected File Format:** Binary RDS symmetric matrix.
- **Header & Structure:** N x N matrix containing pedigree kinship coefficients.
- **Value Types & Encodings:** Numeric covariance.
- **How the Code Processes It:** Covariance structure K in pedigree-based heritability estimation.

### 3. `M`
- **Data Description:** Genotype dosage matrix for genomic comparison.
- **Expected File Format:** Numeric matrix.
- **Header & Structure:** N individuals x M SNPs.
- **Value Types & Encodings:** 0, 1, 2 dosages.
- **How the Code Processes It:** Used to calculate genomic covariance matrix G for side-by-side heritability comparison.

### 4. `pedigree.txt`
- **Data Description:** Raw pedigree records (ID, Sire, Dam).
- **Expected File Format:** Tab-delimited text file.
- **Header & Structure:** 3 columns with header.
- **Value Types & Encodings:** Character IDs.
- **How the Code Processes It:** Fallback source for pedigree matrix calculation if A_matrix.rds is absent.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/pblup_heritability_<trait>.rds: PBLUP variance components and h^2 estimates.`
- `results/gblup_pblup_comparison_<trait>.rds: Comparative summary object.`
- `results/gblup_pblup_heritability_summary.csv: Tabular comparison of pedigree vs genomic heritability.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `rrBLUP`, `AGHmatrix`

To install any missing packages in R:
```R
required_packages <- c("rrBLUP", "AGHmatrix")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 11_pblup_vs_gblup_heritability.R
```
or inside an R interactive session:
```R
source("11_pblup_vs_gblup_heritability.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
