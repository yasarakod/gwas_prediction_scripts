# GWAS & Genomic Prediction Pipeline Repository

This repository provides an organized and fully documented suite of **R scripts** for Genome-Wide Association Study (GWAS) marker prioritization, multi-model genomic prediction (GBLUP, PBLUP, Bayesian regression models, machine learning), heritability estimation, and statistical benchmarking.

---

## 📁 Repository Structure

```
GWAS Prediction/
├── README.md                                    # Comprehensive project documentation
├── 01 Pedigree Matrix Calculation/              # Module 01: Pedigree A & Ainv matrix calculation
│   ├── 01_pedigree_matrix_calculation.R
│   ├── 01_pedigree_matrix_calculation.md
│   └── sample_inputs/                           # Local sample input data
├── 02 Pedigree Data Preprocessing/               # Module 02: Raw data conversion & matrix assembly
│   ├── 02_pedigree_data_preprocessing.R
│   ├── 02_pedigree_data_preprocessing.md
├── 03 Genomic Prediction Cross Validation/      # Module 03: Multi-model CV (GBLUP, Bayes, RF, EN, etc.)
│   ├── 03_genomic_prediction_cross_validation.R
│   ├── 03_genomic_prediction_cross_validation.md
├── 04 Genomic Prediction Top Markers/           # Module 04: 10-fold CV with top GWAS markers
│   ├── 04_genomic_prediction_top_markers.R
│   ├── 04_genomic_prediction_top_markers.md
├── 05 Genomic Prediction CLI Module/            # Module 05: Argparse command-line prediction runner
│   ├── 05_genomic_prediction_cli_module.R
│   ├── 05_genomic_prediction_cli_module.md
├── 06 Weighted GBLUP Prediction/                # Module 06: Iterative Weighted GBLUP (wGBLUP)
│   ├── 06_weighted_gblup_prediction.R
│   ├── 06_weighted_gblup_prediction.md
├── 07 VRMA Prediction Models/                   # Module 07: Final integrated VRMA prediction models
│   ├── 07_vrma_prediction_models.R
│   ├── 07_vrma_prediction_models.md
├── 08 VRMA Prediction Excel Export/             # Module 08: Tabular export of CV iterations
│   ├── 08_vrma_prediction_excel_export.R
│   ├── 08_vrma_prediction_excel_export.md
├── 09 GWAS Markers Population Sizes Prediction/ # Module 09: Training population size benchmarking
│   ├── 09_gwas_markers_population_sizes_prediction.R
│   ├── 09_gwas_markers_population_sizes_prediction.md
├── 10 GBLUP Heritability Estimation/            # Module 10: GBLUP narrow-sense heritability (h2)
│   ├── 10_gblup_heritability_estimation.R
│   ├── 10_gblup_heritability_estimation.md
├── 11 PBLUP vs GBLUP Heritability/              # Module 11: PBLUP vs GBLUP heritability comparison
│   ├── 11_pblup_vs_gblup_heritability.R
│   ├── 11_pblup_vs_gblup_heritability.md
├── 12 Multivariate Genetic Correlation/         # Module 12: Multi-trait genetic/phenotypic correlation
│   ├── 12_multivariate_genetic_correlation.R
│   ├── 12_multivariate_genetic_correlation.md
├── 13 Kaplan Meier Survival Analysis/           # Module 13: Kaplan-Meier survival curves & statistics
│   ├── 13_kaplan_meier_survival_analysis.R
│   ├── 13_kaplan_meier_survival_analysis.md
├── 14 GWAS vs Random Marker Comparison/         # Module 14: GWAS vs random markers (Friedman tests)
│   ├── 14_gwas_vs_random_marker_comparison.R
│   ├── 14_gwas_vs_random_marker_comparison.md
├── 15 Population Size Friedman Analysis/        # Module 15: Population size Friedman testing & plots
│   ├── 15_population_size_friedman_analysis.R
│   ├── 15_population_size_friedman_analysis.md
├── 16 Population Size Statistical Analysis/     # Module 16: Advanced sample size stats (alpha = 0.01)
│   ├── 16_population_size_statistical_analysis.R
│   ├── 16_population_size_statistical_analysis.md
├── 17 Population Size Fixed YAxis Plots/        # Module 17: Standardized y-axis boxplot formatting
│   ├── 17_population_size_fixed_yaxis_plots.R
│   ├── 17_population_size_fixed_yaxis_plots.md
├── 18 GBLUP PopSize Folds Trajectories/         # Module 18: Mean performance trajectory line plots
│   ├── 18_gblup_popsize_folds_trajectories.R
│   ├── 18_gblup_popsize_folds_trajectories.md
├── 19 Summary Graphs Methods Markers PopSizes/  # Module 19: Multi-factor summary visualizations
│   ├── 19_summary_graphs_methods_markers_popsizes.R
│   ├── 19_summary_graphs_methods_markers_popsizes.md
├── 20 Friedman Nemenyi Analysis/                # Module 20: Non-parametric Friedman-Nemenyi tests
│   ├── 20_friedman_nemenyi_analysis.R
│   ├── 20_friedman_nemenyi_analysis.md
├── 21 Enhanced Friedman Ranking Analysis/       # Module 21: Enhanced Friedman CLD analysis & reports
│   ├── 21_enhanced_friedman_ranking_analysis.R
│   ├── 21_enhanced_friedman_ranking_analysis.md
├── 22 VRMA GWAS PBLUP Bar Graphs/               # Module 22: Grouped bar charts with Kendall W
│   ├── 22_vrma_gwas_pblup_bar_graphs.R
│   ├── 22_vrma_gwas_pblup_bar_graphs.md
├── 23 VRMA GWAS PBLUP Bar Graphs Final/         # Module 23: Polished publication bar charts (Viridis)
│   ├── 23_vrma_gwas_pblup_bar_graphs_final.R
│   ├── 23_vrma_gwas_pblup_bar_graphs_final.md
├── 24 VRMA GWAS PBLUP Comprehensive Bar Graphs/ # Module 24: Integrated prediction + bar graph pipeline
│   ├── 24_vrma_gwas_pblup_comprehensive_bar_graphs.R
│   ├── 24_vrma_gwas_pblup_comprehensive_bar_graphs.md
├── 25 Traits Mean SD Bar Graphs/                # Module 25: Mean +- SD bar graphs for 1000 Mb markers
│   ├── 25_traits_mean_sd_bar_graphs.R
│   ├── 25_traits_mean_sd_bar_graphs.md
├── 26 Within Method Trait Comparison/           # Module 26: Within-model trait comparison heatmaps
│   ├── 26_within_method_trait_comparison.R
│   ├── 26_within_method_trait_comparison.md
├── 27 Revised Publication Grouped Bar Plots/    # Module 27: Publication grouped bar charts
│   ├── 27_revised_publication_grouped_bar_plots.R
│   ├── 27_revised_publication_grouped_bar_plots.md
└── 28 Dual Series Random vs GWAS Boxplot/       # Module 28: Dual-series Random vs GWAS boxplots
    ├── 28_dual_series_random_vs_gwas_boxplot.R
    ├── 28_dual_series_random_vs_gwas_boxplot.md
```