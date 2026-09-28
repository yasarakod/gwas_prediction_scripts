# ============================================================================== #
# Script: 13_kaplan_meier_survival_analysis.R
# Title: Kaplan-Meier Survival Analysis for Disease Challenge
# Description:
#   Performs non-parametric Kaplan-Meier survival curve estimation, computes median survival times and survival statistics, and generates survival curve plots.
#
# Inputs:
#   - data/pheno (Phenotype data containing survival time and censoring indicators)
#
# Outputs:
#   - results/kaplan_meier_dpc_time.png (Survival probability curve over hours/time)
#   - results/kaplan_meier_dpc_date.png (Survival probability curve over days)
#
# Required R Packages: survival, survminer, dplyr
# ============================================================================== #

# ------------------------------------------------------------------------------
# Directory Setup & Output Path Configuration
# ------------------------------------------------------------------------------
# Automatically create "results" folder for output deliverables
results_dir <- file.path(getwd(), "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

# Automatically locate shared input data directory
if (dir.exists(file.path(getwd(), "data"))) {
  data_dir <- file.path(getwd(), "data")
} else if (dir.exists(file.path(dirname(getwd()), "data"))) {
  data_dir <- file.path(dirname(getwd()), "data")
} else {
  data_dir <- getwd()
}
data_path <- paste0(data_dir, "/")
# ------------------------------------------------------------------------------


# Load required libraries
library(survival)
library(survminer)
library(dplyr)

# Read the data
pheno <- readRDS("pheno")

# Create survival objects
# For DPC_time (continuous time to event)
surv_obj_time <- Surv(time = pheno$DPC_time, event = pheno$Is_Dead)

# For DPC_date (discrete days to event)
surv_obj_date <- Surv(time = pheno$DPC_date, event = pheno$Is_Dead)

# Fit Kaplan-Meier survival curves
km_fit_time <- survfit(surv_obj_time ~ 1, data = pheno)
km_fit_date <- survfit(surv_obj_date ~ 1, data = pheno)

# Create Kaplan-Meier plots
# Plot for DPC_time
p1 <- ggsurvplot(
  km_fit_time,
  data = pheno,
  title = "Kaplan-Meier Survival Curve - DPC Time",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ggtheme = theme_minimal()
)

# Plot for DPC_date
p2 <- ggsurvplot(
  km_fit_date,
  data = pheno,
  title = "Kaplan-Meier Survival Curve - DPC Date",
  xlab = "Days Post Challenge",
  ylab = "Survival Probability",
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ggtheme = theme_minimal()
)

# Display plots
print(p1)
print(p2)

# Optional: Save plots
ggsave("kaplan_meier_dpc_time.png", plot = p1$plot, width = 10, height = 8, dpi = 300)
ggsave("kaplan_meier_dpc_date.png", plot = p2$plot, width = 10, height = 8, dpi = 300)

# Print summary statistics
cat("Summary for DPC_time:\n")
print(summary(km_fit_time))
cat("\nSummary for DPC_date:\n")
print(summary(km_fit_date))
