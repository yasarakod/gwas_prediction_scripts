# GBLUP Genomic Heritability Estimation

**Script File:** `10_gblup_heritability_estimation.R`  
**Location:** `10 GBLUP Heritability Estimation/10_gblup_heritability_estimation.R`

---

## 1. Brief Description
Estimates genomic narrow-sense heritability (h^2 = sigma_a^2 / (sigma_a^2 + sigma_e^2)) for mortality and morphological traits using genomic relationship matrices via restricted maximum likelihood (REML).

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `02_pedigree_data_preprocessing.R`
- **Module Folder:** `02 Pedigree Data Preprocessing`
- **Role in Pipeline:** Generates M, pheno, and selected.txt.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `M`
- **Data Description:** Genotype dosage matrix (N individuals x M SNPs).
- **Expected File Format:** Numeric matrix object.
- **Header & Structure:** Rows = individual IDs, Columns = SNP marker IDs.
- **Value Types & Encodings:** 0, 1, 2 dosages.
- **How the Code Processes It:** Used to compute VanRaden genomic relationship matrix: G = (Z %*% t(Z)) / (2 * sum(p * (1-p))).

### 2. `pheno / phenotypes.csv`
- **Data Description:** Phenotypic records across mortality and biometric traits.
- **Expected File Format:** Dataframe or CSV file.
- **Header & Structure:** Columns: Is_Dead, DPC_time, DPC_date, Length, Weight, Width.
- **Value Types & Encodings:** Continuous and binary observations.
- **How the Code Processes It:** Response vector y in rrBLUP::mixed.solve(y = y, K = G) to partition additive genetic variance (Va) and residual environmental variance (Ve).

### 3. `selected.txt`
- **Data Description:** List of vetted sample IDs to retain in analysis.
- **Expected File Format:** Plain text table (FID, IID).
- **Header & Structure:** Space-delimited identifiers.
- **Value Types & Encodings:** Character strings.
- **How the Code Processes It:** Filters phenotype and genotype records.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/heritability_<trait>.rds: Variance components (Va, Ve) and narrow-sense heritability point estimates.`
- `results/heritability_summary.csv: Summary table of heritability estimates across all evaluated traits.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `rrBLUP`, `data.table`, `SNPRelate`, `snpReady`

To install any missing packages in R:
```R
required_packages <- c("rrBLUP", "data.table", "SNPRelate", "snpReady")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 10_gblup_heritability_estimation.R
```
or inside an R interactive session:
```R
source("10_gblup_heritability_estimation.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
