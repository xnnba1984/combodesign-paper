# Reproduction Materials
This repository contains the analysis and simulation code that supports all empirical results in the manuscript "Synergy-Informed Design of Platform Trials for Combination Therapies: False Positive Control, Allocation Optimization, and Sample Size Determination". Each script is independent of the ```combodesign``` R package and no external functions are required.

File | Purpose | Key outputs
-----|-----|-----|
multiplicity_control.R | Simulates two correlated test statistics under the global null and compares different multiplicity control methods: NoAdj, Bonferroni, Holm, Dunnett for three error metrics (FWER, FMER, MSFP). | Table summary of false-positive rates; Figure 5. ​
power_analysis.R | Optimizes allocation, derives critical and p values for FWER/FM​ER/MSFP, and finds the minimal total N that yields ≥80 % power across grids of synergy (s) and endpoint correlation (ρ). | Figure 7A (sample-size curves); Figure 7B (optimal allocation). ​
gen_dunnett.R | Calculates critical values by generalized Dunnett's procedure to control multiple testing; computes single-test p-value thresholds for the three error criteria. | Figure 6. ​
PDX_analysis.R | Analyzes the PDX dataset: extracts arm-level correlations, estimates normalized effect sizes & synergy, and applies the allocation optimization to each drug pair. | Table 1 (false-positive control) and Table 2 (optimal N & allocation). ​
false_positive.R | Investigates how arm correlations and allocation ratios shape FWER/FM​ER/MSFP when no multiplicity adjustment is applied. | Figures 3 and 4. ​

## Prerequisites
R ≥ 4.1

R packages: MASS, mvtnorm, Matrix, dplyr, tidyr, ggplot2, patchwork, ggrepel, compiler, ggh4x, readxl, purrr

Install them with
```r
install.packages(c(
  "MASS", "mvtnorm", "Matrix", "dplyr", "tidyr",
  "ggplot2", "patchwork", "ggrepel", "compiler",
  "ggh4x", "readxl", "purrr"
))
```

## Contact
For any questions or issue reports, pleasae open an issue or email nxi@ucla.edu.
