# combodesign-paper

This repository contains the analysis and simulation code that supports all empirical results in the manuscript "Synergy-Informed Design of Platform Trials for Combination Therapies: False Positive Control, Allocation Optimization, and Sample Size Determination". Each script is independent of the ```combodesign``` R package—no external functions are required.

File | Purpose | Key outputs
-----|-----|-----|
multiplicity_control.R | Simulates two correlated test statistics under the global null and compares NoAdj, Bonferroni, Holm, Dunnett for three error metrics (FWER, FMER, MSFP). | CSV summary of false-positive rates; three-panel plot vs correlation. ​
power_analysis.R | Optimises allocation, derives critical values for FWER/FM​ER/MSFP, and finds the minimal total N that yields ≥80 % power across grids of synergy (s) and endpoint correlation (ρ). | Figure 2 (sample-size curves); Figure 3 (optimal allocation). ​
gen_dunnett.R | Generalises the Dunnett critical value to correlated two-arm contrasts; computes single-test p-value thresholds for the three error criteria. | Helper functions find_cstar_*, tables of c\*c^\*c\* and thresholds. ​
PDX_analysis.R | Re-analyses the 41591_2015_BFnm3954 patient-derived xenograft data set: extracts arm-level correlations, estimates normalised effect sizes & synergy, and applies the optimisation routine to each drug pair. | Table 1 (false-positive probabilities) and Table 2 (optimal N & allocation) saved to result/. ​
false_positive.R | Investigates how allocation ratios (1:1:1, 2:1:1, …) and arm correlations shape FWER/FM​ER/MSFP when no multiplicity adjustment is applied. | Heat-map style plots of error inflation across correlation grids. ​

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
For any questions or issue report, pleasae open an issue or e-mail nxi@ucla.edu.
