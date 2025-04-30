###############################################
## 0) Load Packages
###############################################
library(MASS)      # For mvrnorm()
library(dplyr)     # Data handling
library(tidyr)
library(ggplot2)
library(ggrepel)

###############################################
## Check 2x2 correlation validity
###############################################
# We have arms: A, B, AB
# Suppose we have cor(A,B) = rhoA_B, cor(AB,A)=rhoAB_A, cor(AB,B)=rhoAB_B.
# We form the 3x3 correlation matrix among the arms:
#    [ 1        rhoA_B   rhoAB_A ]
#    [ rhoA_B   1        rhoAB_B ]
#    [ rhoAB_A  rhoAB_B  1       ]
# Then check if it's positive semidefinite (PSD).
check_correlation_matrix_valid <- function(rhoAB_A, rhoA_B, rhoAB_B) {
  Cmat <- matrix(c(
    1,       rhoA_B,   rhoAB_A,
    rhoA_B,  1,         rhoAB_B,
    rhoAB_A, rhoAB_B,   1
  ), nrow=3, byrow=TRUE)
  
  # Eigenvalues must be >= 0 for PSD
  ev <- eigen(Cmat, symmetric=TRUE)$values
  all(ev > -1e-12)  # a small tolerance for numerical issues
}

###############################################
## Compute test-statistic correlation rho
###############################################
calc_rho12 <- function(nA, nB, nAB, rhoAB_A, rhoAB_B, rhoA_B, sigma2 = 1) {
  VA  <- sigma2 / nA
  VB  <- sigma2 / nB
  VAB <- sigma2 / nAB
  
  # Covariances
  CAB_A <- rhoAB_A * sigma2 / sqrt(nAB*nA)
  CAB_B <- rhoAB_B * sigma2 / sqrt(nAB*nB)
  
  # Numerator
  num <- CAB_B - CAB_A + VA
  # Denominator
  den <- sqrt( (VAB + VA - 2*CAB_A) * (VB + VA) )
  
  rho <- num / den
  return(rho)
}

###############################################
## 3) Single-scenario simulation for FWER, FMER, MSFP
###############################################
simulate_scenario <- function(nA, nB, nAB, 
                              rhoAB_A, rhoAB_B, rhoA_B,
                              B = 5000, alpha = 0.05) {
  # (A) Check 3x3 correlation among arms
  # is_valid <- check_correlation_matrix_valid(rhoAB_A, rhoA_B, rhoAB_B)
  # if (!is_valid) {
  #   return(data.frame(
  #     nA=nA, nB=nB, nAB=nAB,
  #     rhoAB_A=rhoAB_A, rhoA_B=rhoA_B, rhoAB_B=rhoAB_B,
  #     valid=FALSE, FWER=NA, FMER=NA, MSFP=NA
  #   ))
  # }
  
  # (B) Compute correlation of T1, T2
  r12 <- calc_rho12(nA, nB, nAB, rhoAB_A, rhoAB_B, rhoA_B)
  
  # (C) Check if r12 is NA or out of [-1,1]
  if (is.na(r12) || is.nan(r12) || abs(r12) > 1) {
    return(data.frame(
      nA=nA, nB=nB, nAB=nAB,
      rhoAB_A=rhoAB_A, rhoA_B=rhoA_B, rhoAB_B=rhoAB_B,
      valid=FALSE, FWER=NA, FMER=NA, MSFP=NA
    ))
  }
  
  # (D) Simulate bivariate normal for (T1,T2)
  Sigma <- matrix(c(1, r12, r12, 1), 2, 2)
  Tvals <- MASS::mvrnorm(n=B, mu=c(0,0), Sigma=Sigma)
  T1 <- Tvals[,1]
  T2 <- Tvals[,2]
  
  # (E) Rejection logic
  zcrit <- qnorm(1 - alpha/2)
  reject1 <- abs(T1) > zcrit
  reject2 <- abs(T2) > zcrit
  
  # (F) Compute FWER, FMER, MSFP
  fwer <- mean(reject1 | reject2)
  fmer <- mean(reject1 & reject2)
  msfp <- mean((T1 > zcrit) & (T2 > zcrit))
  
  # Return valid scenario
  data.frame(
    nA=nA, nB=nB, nAB=nAB,
    rhoAB_A=rhoAB_A, rhoA_B=rhoA_B, rhoAB_B=rhoAB_B,
    valid=TRUE,
    FWER=fwer, FMER=fmer, MSFP=msfp
  )
}


###############################################
##    Vary one correlation at a time
##    across 0..1 step 0.01, fix the other two at baseline
##    Do for each of the 4 allocation ratio patterns
###############################################

# base sample size multiplier
k <- 100

# Define the 4 ratio patterns
alloc_scenarios <- list(
  c(1,1,1),
  c(2,1,1),
  c(1,2,1),
  c(1,1,2)
)
# We'll fix the other two correlations at 0.3 as an example baseline
rho_baseline <- 0.3

# We'll create a grid for each correlation parameter, one at a time.
rho_grid <- seq(0, 1, by=0.01)

# Prepare a container for final results
all_results <- list()

# Vary rhoAB_A, fix rhoAB_B at 0.3
res_rhoAB_A <- list()
for (alloc in alloc_scenarios) {
  nA   <- alloc[1]*k
  nB   <- alloc[2]*k
  nAB  <- alloc[3]*k
  
  tmp_list <- lapply(rho_grid, function(rval) {
    simulate_scenario(
      nA=nA, nB=nB, nAB=nAB,
      rhoAB_A=rval,       # varied
      rhoA_B=0.01, 
      rhoAB_B=rho_baseline,
      B=100000, alpha=0.05
    )
  })
  df_out <- do.call(rbind, tmp_list)
  # add columns to label what factor we varied
  df_out$varied_factor  <- "rhoAB_A"
  df_out$allocation_str <- paste(alloc, collapse=":")
  res_rhoAB_A[[ paste(alloc, collapse="_") ]] <- df_out
}
res_rhoAB_A <- do.call(rbind, res_rhoAB_A)
all_results[["rhoAB_A"]] <- res_rhoAB_A

# Vary rhoAB_B, fix rhoAB_A and rhoA_B at 0.3
res_rhoAB_B <- list()
for (alloc in alloc_scenarios) {
  nA   <- alloc[1]*k
  nB   <- alloc[2]*k
  nAB  <- alloc[3]*k
  
  tmp_list <- lapply(rho_grid, function(rval) {
    simulate_scenario(
      nA=nA, nB=nB, nAB=nAB,
      rhoAB_A=rho_baseline,
      rhoA_B=0.01,
      rhoAB_B=rval,        # varied
      B=100000, alpha=0.05
    )
  })
  df_out <- do.call(rbind, tmp_list)
  df_out$varied_factor  <- "rhoAB_B"
  df_out$allocation_str <- paste(alloc, collapse=":")
  res_rhoAB_B[[ paste(alloc, collapse="_") ]] <- df_out
}
res_rhoAB_B <- do.call(rbind, res_rhoAB_B)
all_results[["rhoAB_B"]] <- res_rhoAB_B

###############################################
## 5) Combine Everything Into One Data Frame
###############################################
final_df <- do.call(rbind, all_results)

# Drop the invalid scenarios
final_df_valid <- final_df %>% filter(valid)

###############################################
##    For each factor varied, and each allocation ratio,
##    plot FWER vs. the varied correlation
###############################################

plot_df <- final_df_valid %>%
  filter(varied_factor == "rhoAB_B", valid == TRUE)

plot_df_long <- plot_df %>%
  pivot_longer(
    cols = c("FWER", "FMER", "MSFP"),
    names_to = "Metric",
    values_to = "Rate"
  )

plot_df_long$allocation_str <- factor(
  plot_df_long$allocation_str,
  levels = c("1:1:1","2:1:1","1:2:1","1:1:2")
)

plot_df_long$Metric <- factor(
  plot_df_long$Metric,
  levels = c("FWER","FMER","MSFP")
)

# Step 2a) Identify all unique row/column facet values:
unique_allocs <- unique(plot_df_long$allocation_str)
unique_metrics <- unique(plot_df_long$Metric)

# Step 2b) Expand to get every combination
baseline_df <- expand.grid(
  allocation_str = unique_allocs,
  Metric         = unique_metrics
)

# Step 2c) Assign the baseline "independent" error rate per metric
# For the standard 2-test scenario with alpha=0.05:
baseline_df$baseline_value <- NA_real_
baseline_df$baseline_value[baseline_df$Metric=="FWER"] <- 0.0975
baseline_df$baseline_value[baseline_df$Metric=="FMER"] <- 0.0025
baseline_df$baseline_value[baseline_df$Metric=="MSFP"] <- 0.000625

df_xloc <- plot_df_long %>%
  group_by(allocation_str, Metric) %>%
  summarise(
    x_loc = max(rhoAB_B, na.rm = TRUE),
    .groups = "drop"
  )

# Merge that into baseline_df
df_labels <- baseline_df %>%
  left_join(df_xloc, by = c("allocation_str","Metric")) %>%
  # Format the label: "Baseline=0.0975" etc.
  mutate(label_txt = paste0("Baseline=", round(baseline_value, 6)))


ggplot(plot_df_long, aes(x = rhoAB_B, y = Rate)) +
  geom_line(size = 1.05, color = "black") +
  ggh4x::facet_grid2(
    rows = vars(allocation_str),
    cols = vars(Metric),
    scales = "free_y", independent = "y"
  ) +
  # White background, bigger base font, etc.
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10) 
  ) +
  labs(
    x = expression(rho[AB*","*B]),
    y = "False Positive Rate"
  ) +
  # Add a second axis on the right side, labeled in math format
  scale_y_continuous(
    sec.axis = sec_axis(
      transform  = ~ .,   # identity transform
      name   = expression(n[A]~":"~n[B]~":"~n[AB]),
      breaks = NULL,  # hide numeric ticks
      labels = NULL
    )
  ) +
  geom_hline(
    data        = baseline_df,
    aes(yintercept = baseline_value),
    color       = "blue",
    linetype    = "dashed",
    size        = 0.6,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))











