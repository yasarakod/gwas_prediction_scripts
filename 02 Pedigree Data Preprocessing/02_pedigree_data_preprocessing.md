# Pedigree & Genomic Data Preprocessing Pipeline

**Script File:** `02_pedigree_data_preprocessing.R`  
**Location:** `02 Pedigree Data Preprocessing/02_pedigree_data_preprocessing.R`

---

## 1. Brief Description
End-to-end preprocessing workflow that cleans raw pedigree files, converts PLINK raw/ped/map genotypes into numerical dosage matrices, aligns phenotypic records, and prepares analysis-ready R objects.

---

## 2. Pre-required Scripts (Upstream Pipeline Dependencies)
To generate the exact input files required for this script, the following upstream scripts should be executed in order (or verify their presence in `sample_inputs/` or `data/`):

### 1. `None (Foundational Preprocessing Script)`
- **Role in Pipeline:** Operates directly on raw PLINK binary/text outputs (tempgeno.raw, tempgeno.map, tempgeno.ped), raw pedigree (pedigree.txt), and phenotype records (phenotypes.csv). No prior R scripts are required.

> **Execution Shortcut:** If pre-generated datasets or previous run outputs are already present in the local `sample_inputs/` directory or `../data/`, you do **not** need to re-run the upstream scripts and can execute this script immediately.

---

## 3. Inputs (Detailed Specification)
The script expects the following input files. For each file, the required biological/statistical data, expected file format, structure, and code processing expectations are detailed below:

### 1. `pedigree.txt`
- **Data Description:** Individual parentage table (ID, Sire, Dam) used to construct ancestral pedigree depth.
- **Expected File Format:** Tab-delimited text file with header.
- **Header & Structure:** Columns: ID, Sire, Dam.
- **Value Types & Encodings:** Alphanumeric sample IDs; 0 or NA for unknown parents.
- **How the Code Processes It:** Processed to impute missing parents, verify topological ordering, and format into completed_pedigree.rds.

### 2. `phenotypes.csv`
- **Data Description:** Experimental challenge phenotype records for all evaluated individuals.
- **Expected File Format:** Comma-separated values (CSV) with header.
- **Header & Structure:** Columns: ID, Is_Dead, DPC_time, DPC_date, Length, Weight, Width.
- **Value Types & Encodings:** ID: Character string matching genotype IDs; Is_Dead: Binary integer (1=Dead, 0=Alive); DPC_time: Numeric hours to mortality; DPC_date: Numeric days to mortality; Biometrics: Floating-point measurements.
- **How the Code Processes It:** Read via read.csv(), filtered to match genotyped individuals, cleaned of outliers, and serialized to pheno.

### 3. `tempgeno.raw`
- **Data Description:** PLINK additive genotype dosage table generated from plink --recodeA.
- **Expected File Format:** Space-delimited text file with header.
- **Header & Structure:** Leading columns: FID, IID, PAT, MAT, SEX, PHENOTYPE; followed by SNP columns named SNP_Allele (e.g. AX-123456_A).
- **Value Types & Encodings:** Additive integer dosages: 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate, NA = missing.
- **How the Code Processes It:** Extracted into dosage matrix M (rows = individuals, columns = SNPs) via data.table::fread().

### 4. `tempgeno.map`
- **Data Description:** PLINK variant map file detailing genomic positions.
- **Expected File Format:** Tab-delimited 4-column map file without header.
- **Header & Structure:** Column 1: Chromosome / Contig; Column 2: Variant ID (rs / SNP name); Column 3: Genetic distance (cM, 0 if unmapped); Column 4: Physical base-pair position (bp).
- **Value Types & Encodings:** Chromosomes as integer/string; positions as positive integer base pairs.
- **How the Code Processes It:** Loaded to assemble hapmap dataframe for marker quality control and chromosome alignment.

### 5. `selected.txt / selected.fam`
- **Data Description:** Whitelist of sample identifiers meeting sample-level call rate and pedigree inclusion criteria.
- **Expected File Format:** Space-delimited plain text file without header.
- **Header & Structure:** Column 1: Family ID (FID); Column 2: Individual ID (IID).
- **Value Types & Encodings:** Exact character identifiers matching PLINK .fam files.
- **How the Code Processes It:** Subsets raw genotype matrices to retain only vetted high-quality samples.

> **Note on Data Locations:** Sample input files for quick execution are bundled directly in the local `sample_inputs/` folder within this directory. Global datasets are also available in the top-level `../data/` folder.

---

## 4. Outputs
The script automatically creates and writes its deliverables into the dedicated `results/` folder:
- `results/M: Numeric genotype dosage matrix (individuals x SNPs).`
- `results/pheno: Formatted phenotype dataframe matched to genotype rows.`
- `results/hapmap: Variant annotation dataframe with marker names, chromosomes, and physical coordinates.`
- `results/A_matrix.rds: Numerator additive relationship matrix computed from pedigree.`

---

## 5. Dependencies & Required Packages
Please ensure the following R packages are installed before execution:
- `dplyr`, `data.table`, `AGHmatrix`, `SNPRelate`, `snpReady`

To install any missing packages in R:
```R
required_packages <- c("dplyr", "data.table", "AGHmatrix", "SNPRelate", "snpReady")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
```

---

## 6. Execution Instructions
Run the script from within its folder or from the project root directory:
```bash
# From within the module directory:
Rscript 02_pedigree_data_preprocessing.R
```
or inside an R interactive session:
```R
source("02_pedigree_data_preprocessing.R")
```
All results, models, and figures will be saved automatically into the `results/` folder.
