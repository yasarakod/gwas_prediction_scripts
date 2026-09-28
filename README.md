# GWAS & Genomic Prediction Pipeline Repository

This repository provides an organized, deduplicated, and fully documented suite of **R scripts** and **bundled datasets** for Genome-Wide Association Study (GWAS) marker prioritization, multi-model genomic prediction (GBLUP, PBLUP, Bayesian regression models, machine learning), heritability estimation, and statistical benchmarking.

---

## 📁 Repository Structure

```
GWAS Prediction/
├── README.md                                    # Comprehensive project documentation
├── data/                                        # Consolidated input datasets (genotypes, phenotypes, pedigree, etc.)
│   ├── phenotypes.csv                           # Phenotypic observations
│   ├── pedigree.txt                             # Pedigree records (ID, Sire, Dam)
│   ├── genotype_data.bed/.bim/.fam              # PLINK binary genotype files
│   ├── wald_*.txt                               # GWAS Wald test summary statistics
│   ├── top_*_snps.txt                           # Prioritized marker subsets
│   └── ...                                      # Preprocessed relationship matrices & dosage data
│
├── 01 Pedigree Matrix Calculation/              # Module 01: Pedigree A & Ainv matrix calculation
│   ├── 01_pedigree_matrix_calculation.R
│   ├── 01_pedigree_matrix_calculation.md
│   └── sample_inputs/                           # Local sample input data
├── 02 Pedigree Data Preprocessing/               # Module 02: Raw data conversion & matrix assembly
│   ├── 02_pedigree_data_preprocessing.R
│   ├── 02_pedigree_data_preprocessing.md
│   └── sample_inputs/
├── 03 Genomic Prediction Cross Validation/      # Module 03: Multi-model CV (GBLUP, Bayes, RF, EN, etc.)
│   ├── 03_genomic_prediction_cross_validation.R
│   ├── 03_genomic_prediction_cross_validation.md
│   └── sample_inputs/
├── 04 Genomic Prediction Top Markers/           # Module 04: 10-fold CV with top GWAS markers
│   ├── 04_genomic_prediction_top_markers.R
│   ├── 04_genomic_prediction_top_markers.md
│   └── sample_inputs/
├── 05 Genomic Prediction CLI Module/            # Module 05: Argparse command-line prediction runner
│   ├── 05_genomic_prediction_cli_module.R
│   ├── 05_genomic_prediction_cli_module.md
│   └── sample_inputs/
├── 06 Weighted GBLUP Prediction/                # Module 06: Iterative Weighted GBLUP (wGBLUP)
│   ├── 06_weighted_gblup_prediction.R
│   ├── 06_weighted_gblup_prediction.md
│   └── sample_inputs/
├── 07 VRMA Prediction Models/                   # Module 07: Final integrated VRMA prediction models
│   ├── 07_vrma_prediction_models.R
│   ├── 07_vrma_prediction_models.md
│   └── sample_inputs/
├── 08 VRMA Prediction Excel Export/             # Module 08: Tabular export of CV iterations
│   ├── 08_vrma_prediction_excel_export.R
│   ├── 08_vrma_prediction_excel_export.md
│   └── sample_inputs/
├── 09 GWAS Markers Population Sizes Prediction/ # Module 09: Training population size benchmarking
│   ├── 09_gwas_markers_population_sizes_prediction.R
│   ├── 09_gwas_markers_population_sizes_prediction.md
│   └── sample_inputs/
├── 10 GBLUP Heritability Estimation/            # Module 10: GBLUP narrow-sense heritability (h2)
│   ├── 10_gblup_heritability_estimation.R
│   ├── 10_gblup_heritability_estimation.md
│   └── sample_inputs/
├── 11 PBLUP vs GBLUP Heritability/              # Module 11: PBLUP vs GBLUP heritability comparison
│   ├── 11_pblup_vs_gblup_heritability.R
│   ├── 11_pblup_vs_gblup_heritability.md
│   └── sample_inputs/
├── 12 Multivariate Genetic Correlation/         # Module 12: Multi-trait genetic/phenotypic correlation
│   ├── 12_multivariate_genetic_correlation.R
│   ├── 12_multivariate_genetic_correlation.md
│   └── sample_inputs/
├── 13 Kaplan Meier Survival Analysis/           # Module 13: Kaplan-Meier survival curves & statistics
│   ├── 13_kaplan_meier_survival_analysis.R
│   ├── 13_kaplan_meier_survival_analysis.md
│   └── sample_inputs/
├── 14 GWAS vs Random Marker Comparison/         # Module 14: GWAS vs random markers (Friedman tests)
│   ├── 14_gwas_vs_random_marker_comparison.R
│   ├── 14_gwas_vs_random_marker_comparison.md
│   └── sample_inputs/
├── 15 Population Size Friedman Analysis/        # Module 15: Population size Friedman testing & plots
│   ├── 15_population_size_friedman_analysis.R
│   ├── 15_population_size_friedman_analysis.md
│   └── sample_inputs/
├── 16 Population Size Statistical Analysis/     # Module 16: Advanced sample size stats (alpha = 0.01)
│   ├── 16_population_size_statistical_analysis.R
│   ├── 16_population_size_statistical_analysis.md
│   └── sample_inputs/
├── 17 Population Size Fixed YAxis Plots/        # Module 17: Standardized y-axis boxplot formatting
│   ├── 17_population_size_fixed_yaxis_plots.R
│   ├── 17_population_size_fixed_yaxis_plots.md
│   └── sample_inputs/
├── 18 GBLUP PopSize Folds Trajectories/         # Module 18: Mean performance trajectory line plots
│   ├── 18_gblup_popsize_folds_trajectories.R
│   ├── 18_gblup_popsize_folds_trajectories.md
│   └── sample_inputs/
├── 19 Summary Graphs Methods Markers PopSizes/  # Module 19: Multi-factor summary visualizations
│   ├── 19_summary_graphs_methods_markers_popsizes.R
│   ├── 19_summary_graphs_methods_markers_popsizes.md
│   └── sample_inputs/
├── 20 Friedman Nemenyi Analysis/                # Module 20: Non-parametric Friedman-Nemenyi tests
│   ├── 20_friedman_nemenyi_analysis.R
│   ├── 20_friedman_nemenyi_analysis.md
│   └── sample_inputs/
├── 21 Enhanced Friedman Ranking Analysis/       # Module 21: Enhanced Friedman CLD analysis & reports
│   ├── 21_enhanced_friedman_ranking_analysis.R
│   ├── 21_enhanced_friedman_ranking_analysis.md
│   └── sample_inputs/
├── 22 VRMA GWAS PBLUP Bar Graphs/               # Module 22: Grouped bar charts with Kendall W
│   ├── 22_vrma_gwas_pblup_bar_graphs.R
│   ├── 22_vrma_gwas_pblup_bar_graphs.md
│   └── sample_inputs/
├── 23 VRMA GWAS PBLUP Bar Graphs Final/         # Module 23: Polished publication bar charts (Viridis)
│   ├── 23_vrma_gwas_pblup_bar_graphs_final.R
│   ├── 23_vrma_gwas_pblup_bar_graphs_final.md
│   └── sample_inputs/
├── 24 VRMA GWAS PBLUP Comprehensive Bar Graphs/ # Module 24: Integrated prediction + bar graph pipeline
│   ├── 24_vrma_gwas_pblup_comprehensive_bar_graphs.R
│   ├── 24_vrma_gwas_pblup_comprehensive_bar_graphs.md
│   └── sample_inputs/
├── 25 Traits Mean SD Bar Graphs/                # Module 25: Mean +- SD bar graphs for 1000 Mb markers
│   ├── 25_traits_mean_sd_bar_graphs.R
│   ├── 25_traits_mean_sd_bar_graphs.md
│   └── sample_inputs/
├── 26 Within Method Trait Comparison/           # Module 26: Within-model trait comparison heatmaps
│   ├── 26_within_method_trait_comparison.R
│   ├── 26_within_method_trait_comparison.md
│   └── sample_inputs/
├── 27 Revised Publication Grouped Bar Plots/    # Module 27: Publication grouped bar charts
│   ├── 27_revised_publication_grouped_bar_plots.R
│   ├── 27_revised_publication_grouped_bar_plots.md
│   └── sample_inputs/
└── 28 Dual Series Random vs GWAS Boxplot/       # Module 28: Dual-series Random vs GWAS boxplots
    ├── 28_dual_series_random_vs_gwas_boxplot.R
    ├── 28_dual_series_random_vs_gwas_boxplot.md
    └── sample_inputs/
```

---

## 🚀 Key Improvements & Pipeline Features
1. **Deduplication:** Removed identical files, redundant copies, and scratch test scripts across the source directories.
2. **Clean Naming with Spaces:** Each script is located in its own clean, space-separated, numbered subfolder.
3. **Detailed Input Documentation:** Every script is paired with a dedicated `<ScriptName>.md` README file providing in-depth data descriptions, file formats, columns/dimensions, value types/encodings, and code usage expectations for every single input file.
4. **Self-Contained Sample Inputs:** Each module contains a `sample_inputs/` folder populated with the exact input files needed to run the script.
5. **Automated Output Handling:** Scripts automatically create and write deliverables into a local `results/` folder.