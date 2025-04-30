library(MASS)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(compiler)

set.seed(12345)

################################################################################
## 1) Define the scenario grids and fixed parameters
################################################################################

synergy_vec <- seq(0.7, 1.3, by = 0.1)  # synergy values
rho_vec     <- seq(0.2, 0.8, by = 0.2)    # grid for "rho" used in simulation of Z's
rho_vec <- c(0.1, 0.3, 0.5, 0.7)

deltaB       <- 0.3       # effect size for B vs. A
target_power <- 0.80      # target power (80%)
nrep_opt     <- 10000      # replications for the optimization stage
nrep_final   <- 50000      # replications for the final check

# Error‐rate definitions for 2 tests:
alpha_fwer <- 0.05      # Family-Wise Error Rate (FWER) control
alpha_fmer <- 0.0025    # FMER control
alpha_msfp <- 0.000625  # MSFP control

################################################################################
## 2) Probability that both |Z1| > c and |Z2| > c for bivariate normal
################################################################################
prob_both_above_c <- function(c, rho) {
  # Compute Probability(Z1 > c, Z2 > c) for bivariate normal (mean 0, correlation=rho)
  # and multiply by 4 to account for all four quadrants.
  pmass_quad <- pmvnorm(
    lower = c(c, c),
    upper = c(Inf, Inf),
    mean  = c(0, 0),
    corr  = matrix(c(1, rho, rho, 1), nrow = 2)
  )
  4 * pmass_quad
}
prob_both_above_c <- cmpfun(prob_both_above_c)

################################################################################
## 3) Get critical value c* for two‐sided tests for each error type
################################################################################
crit_cache <- list()  # Global cache to avoid repeated computations

calc_crit_2sided <- function(rho, error_type="FWER") {
  key <- paste0(error_type, "_", round(rho, 4))
  if (!is.null(crit_cache[[key]])) {
    return(crit_cache[[key]])
  }
  
  if (error_type == "FWER") {
    # Want: P(max(|Z1|,|Z2|) <= c) = 1 - alpha_fwer
    f_c <- function(c) {
      p_in <- pmvnorm(
        lower = c(-c, -c),
        upper = c(c, c),
        mean = c(0, 0),
        corr = matrix(c(1, rho, rho, 1), nrow = 2)
      )
      p_in - (1 - alpha_fwer)
    }
  } else if (error_type == "FMER") {
    # Want: P(|Z1| > c & |Z2| > c) = alpha_fmer
    f_c <- function(c) {
      prob_both_above_c(c, rho) - alpha_fmer
    }
  } else if (error_type == "MSFP") {
    # Want: P(Z1 > c & Z2 > c) = alpha_msfp
    f_c <- function(c) {
      prob_both_above_c(c, rho)  - alpha_msfp
    }
  } else {
    stop("Unknown error_type: must be FWER, FMER, or MSFP")
  }
  
  out <- uniroot(f_c, interval = c(0, 5))
  c_star <- out$root
  crit_cache[[key]] <<- c_star
  c_star
}
calc_crit_2sided <- cmpfun(calc_crit_2sided)

################################################################################
## 4) Correlation of (Z1,Z2) under the global null
################################################################################
get_corr_Z1_Z2_null <- function(pA, pB, pAB, rho, sigma2 = 1, N = 100) {
  nA  <- round(N * pA)
  nB  <- round(N * pB)
  nAB <- round(N * pAB)
  if (min(nA, nB, nAB) < 2) return(NA)
  
  varA  <- sigma2 / nA
  varB  <- sigma2 / nB
  varAB <- sigma2 / nAB
  
  covAB_A <- rho * sqrt(varA * varAB)
  #covB_A  <- 0.001 * sqrt(varA * varB)
  covB_A  <- 0
  covB_AB <- rho * sqrt(varB * varAB)
  
  var_diffAB_A <- varAB + varA - 2 * covAB_A
  var_diffB_A  <- varB + varA - 2 * covB_A
  if (var_diffAB_A <= 0 || var_diffB_A <= 0) return(NA)
  
  cov_diff <- covB_AB - covB_A - covAB_A + varA
  rhoZ <- cov_diff / sqrt(var_diffAB_A * var_diffB_A)
  if (!is.finite(rhoZ) || abs(rhoZ) > 1) return(NA)
  rhoZ
}
get_corr_Z1_Z2_null <- cmpfun(get_corr_Z1_Z2_null)

################################################################################
## 5) Simulate draws (Z1, Z2) from the alternative hypothesis
################################################################################
simulate_draws_Z <- function(pA, pB, pAB,
                             muA, muB, muAB,
                             rho, sigma2 = 1, N = 100,
                             nrep = 2000) {
  nA  <- round(N * pA)
  nB  <- round(N * pB)
  nAB <- round(N * pAB)
  if (min(nA, nB, nAB) < 2) {
    return(matrix(NA, nrow = nrep, ncol = 2))
  }
  
  varA  <- sigma2 / nA
  varB  <- sigma2 / nB
  varAB <- sigma2 / nAB
  
  covAB_A <- rho * sqrt(varA * varAB)
  #covB_A  <- 0.001 * sqrt(varA * varB)
  covB_A  <- 0
  covB_AB <- rho * sqrt(varB * varAB)
  
  var_diffAB_A <- varAB + varA - 2 * covAB_A
  var_diffB_A  <- varB + varA - 2 * covB_A
  if (var_diffAB_A <= 0 || var_diffB_A <= 0) {
    return(matrix(NA, nrow = nrep, ncol = 2))
  }
  
  Sigma_means <- matrix(c(
    varA,      covB_A,    covAB_A,
    covB_A,    varB,      covB_AB,
    covAB_A,   covB_AB,   varAB
  ), nrow = 3, byrow = TRUE)
  
  L <- tryCatch(chol(Sigma_means), error = function(e) NULL)
  if (is.null(L)) {
    return(matrix(NA, nrow = nrep, ncol = 2))
  }
  
  Z <- matrix(rnorm(nrep * 3), nrep, 3)
  # Shift by (muA, muB, muAB)
  draws <- matrix(rep(c(muA, muB, muAB), each = nrep),
                  nrow = nrep, ncol = 3) + Z %*% t(L)
  
  diffAB_A <- draws[, 3] - draws[, 1]
  diffB_A  <- draws[, 2] - draws[, 1]
  z1 <- diffAB_A / sqrt(var_diffAB_A)
  z2 <- diffB_A  / sqrt(var_diffB_A)
  
  cbind(z1, z2)
}
simulate_draws_Z <- cmpfun(simulate_draws_Z)

################################################################################
## 6) Simulate power for a given error_type using the chosen design
################################################################################
simulate_power_2sided <- function(N, pA, pB, pAB,
                                  deltaB, synergy,
                                  rho, error_type = "FWER",
                                  sigma2 = 1, nrep = 2000) {
  rhoZ <- get_corr_Z1_Z2_null(pA, pB, pAB, rho, sigma2, N)
  if (is.na(rhoZ)) return(c(power_AB = NA, power_B = NA))
  
  c_star <- calc_crit_2sided(rhoZ, error_type = error_type)
  
  muA  <- 0
  muB  <- deltaB
  muAB <- synergy * deltaB
  
  zmat <- simulate_draws_Z(pA, pB, pAB,
                           muA, muB, muAB,
                           rho, sigma2, N, nrep)
  
  valid <- complete.cases(zmat)
  if (sum(valid) == 0) return(c(power_AB = NA, power_B = NA))
  
  zmat <- zmat[valid, , drop = FALSE]
  power_AB <- mean(abs(zmat[, 1]) > c_star)
  power_B  <- mean(abs(zmat[, 2]) > c_star)
  
  if(error_type == "MSFP"){
    power_AB <- mean(zmat[, 1] > c_star)
    power_B  <- mean(zmat[, 2] > c_star)
  }
  
  c(power_AB = power_AB, power_B = power_B)
}
simulate_power_2sided <- cmpfun(simulate_power_2sided)

################################################################################
## 7) Optimal Allocation via Continuous Optimization (Numerical Method)
################################################################################
optimal_allocation <- function(synergy, deltaB, rho_ABA = 0.3, rho_AB = 0.001) {
  # Objective function: maximize min(W1, W2) where
  #   p_AB = x, p_A = (1 - x)*y, p_B = 1 - p_AB - p_A,
  # and with the feasibility constraint that all probabilities >= 0.1.
  # The formulas below approximate the noncentrality parameters:
  #   D1 = 1/p_AB + 1/p_A - 2*rho_ABA/sqrt(p_AB*p_A)
  #   W1 = (synergy*deltaB)^2 / D1
  #   D2 = 1/p_B  + 1/p_A - 2*rho_AB/sqrt(p_A*p_B)
  #   W2 = (deltaB)^2 / D2
  # We maximize min(W1, W2) by minimizing its negative.
  objFun <- function(par) {
    x <- par[1]
    y <- par[2]
    pAB <- x
    pA  <- (1 - x) * y
    pB  <- 1 - pAB - pA
    # Enforce feasibility: each allocation must be at least 0.1
    if (pAB < 0.01 || pA < 0.01 || pB < 0.01) return(1e6)
    
    D1 <- (1 / pAB) + (1 / pA) - 2 * rho_ABA / sqrt(pAB * pA)
    D2 <- (1 / pB)  + (1 / pA) - 2 * rho_AB   / sqrt(pA * pB)
    W1 <- (synergy * deltaB)^2 / D1
    W2 <- (deltaB)^2 / D2
    return(-min(W1, W2))
  }
  
  res <- optim(
    par    = c(0.3, 0.3),  # initial guess for (x, y)
    fn     = objFun,
    method = "L-BFGS-B",
    lower  = c(0.0001, 0.0001),
    upper  = c(0.9999, 0.9999)
  )
  
  xOpt <- res$par[1]
  yOpt <- res$par[2]
  pAB  <- xOpt
  pA   <- (1 - xOpt) * yOpt
  pB   <- 1 - pAB - pA
  
  # Return the optimal allocation and the achieved objective value.
  return(c(pA = pA, pB = pB, pAB = pAB, max_minW = -res$value))
}

################################################################################
## 8) Use the allocation method then simulate power
################################################################################
evaluate_N <- function(N, synergy, rho,
                       error_type = "FWER",
                       deltaB = 0.3,
                       nrep_opt = nrep_opt) {
  # Get the optimal allocation using the numerical method.
  alloc <- optimal_allocation(synergy, deltaB, rho_ABA = rho, rho_AB = 0.001)
  
  # Compute simulated power using these allocation probabilities.
  pow <- simulate_power_2sided(N, pA = alloc["pA"], pB = alloc["pB"], pAB = alloc["pAB"],
                               deltaB = deltaB, synergy = synergy,
                               rho = rho, error_type = error_type, nrep = nrep_opt)
  # Use the minimum of the two power values (for the two tests).
  power_val <- min(pow["power_AB"], pow["power_B"])
  list(
    power = power_val,
    allocation = alloc
  )
}
evaluate_N <- cmpfun(evaluate_N)

################################################################################
## 9) Binary search: Find minimal N that achieves target power
################################################################################
find_minN_for_80pct_binary <- function(N_start, synergy, rho,
                                       error_type = "FWER",
                                       deltaB = 0.3,
                                       target_power = 0.8,
                                       nrep_opt = nrep_opt,
                                       nrep_final = nrep_final,
                                       tol = -0.005) {
  # Step 1: Find an upper bound for N
  N_low <- N_start
  eval_low <- evaluate_N(N_low, synergy, rho, error_type, deltaB, nrep_opt)
  if (is.na(eval_low$power)) eval_low$power <- 0
  
  if (eval_low$power >= target_power) {
    N_high <- N_low
  } else {
    N_high <- N_low
    repeat {
      N_high <- N_high * 2
      eval_high <- evaluate_N(N_high, synergy, rho, error_type, deltaB, nrep_opt)
      if (!is.na(eval_high$power) && eval_high$power >= target_power) break
      if (N_high > 5000) break
    }
  }
  
  # Step 2: Binary search for minimal N between N_low and N_high
  while ((N_high - N_low) > 1) {
    N_mid <- floor((N_low + N_high) / 2)
    eval_mid <- evaluate_N(N_mid, synergy, rho, error_type, deltaB, nrep_opt)
    if (!is.na(eval_mid$power) && eval_mid$power >= target_power) {
      N_high <- N_mid
    } else {
      N_low <- N_mid
    }
  }
  
  # Step 3: Final check with increased replication
  final_eval <- evaluate_N(N_high, synergy, rho, error_type, deltaB, nrep_final)
  while ((is.na(final_eval$power) || final_eval$power < target_power) && N_high < 5000) {
    N_high <- N_high + 1
    final_eval <- evaluate_N(N_high, synergy, rho, error_type, deltaB, nrep_final)
  }
  
  # Step 4: Try decrementing N if possible
  candidate <- N_high
  repeat {
    eval_candidate <- evaluate_N(candidate - 1, synergy, rho, error_type, deltaB, nrep_final)
    if (is.na(eval_candidate$power) || eval_candidate$power < target_power - tol)
      break
    candidate <- candidate - 1
    if (candidate <= N_start) break
  }
  
  final_evals <- evaluate_N(candidate, synergy, rho, error_type, deltaB, nrep_final)
  
  list(
    N = candidate,
    ratio = final_evals$allocation,
    power = final_evals$power
  )
}
find_minN_for_80pct_binary <- cmpfun(find_minN_for_80pct_binary)


################################################################################
## 10) Main Loop: Loop over error types, synergy, and rho to determine N
################################################################################
params <- expand.grid(
  synergy = synergy_vec,
  rho     = rho_vec
)
results_list <- list()
idx <- 1

#start_time <- Sys.time()

# Loop over error_type choices and each (synergy, rho) combination:
for (error_type in c("FWER", "FMER", "MSFP")) {
  for (i in seq_len(nrow(params))) {
    s_val <- params$synergy[i]
    r_val <- params$rho[i]
    
    cat(sprintf("error_type=%s | synergy=%.2f | rho=%.2f...\n",
                error_type, s_val, r_val))
    
    out <- find_minN_for_80pct_binary(
      N_start = 20,
      synergy = s_val,
      rho = r_val,
      error_type = error_type,
      deltaB = deltaB,
      target_power = target_power,
      nrep_opt = nrep_opt,
      nrep_final = nrep_final,
      tol = -0.005
    )
    
    results_list[[idx]] <- data.frame(
      error_type = error_type,
      synergy    = s_val,
      rho        = r_val,
      N          = out$N,
      pA_opt     = out$ratio["pA"],
      pB_opt     = out$ratio["pB"],
      pAB_opt    = out$ratio["pAB"],
      power      = out$power
    )
    idx <- idx + 1
  }
}

df_results <- do.call(rbind, results_list)
#end_time <- Sys.time()

# cat("\nFinal Results:\n")
# print(df_results)
# cat("\nTotal Running Time:", end_time - start_time, "\n")


library(ggplot2)
library(dplyr)
library(tidyr)

# Define a custom color palette 
custom_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

# Ensure that 'rho' is a factor with levels 0.2, 0.4, 0.6, 0.8 (even if some are missing in the data)
df_results <- df_results %>%
  mutate(rho = factor(rho, levels = c(0.1, 0.3, 0.5, 0.7)))

df_results$error_type <- factor(df_results$error_type, levels = c("FWER", "FMER", "MSFP"))

y_min <- floor(min(df_results$N) / 100) * 100
y_max <- ceiling(max(df_results$N) / 100) * 100

### First Plot: 3x1 Facet Plot (Total Sample Size vs. Synergy) with unified y-axis ###
ggplot(df_results, aes(x = synergy, y = N, color = rho, group = rho)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~ error_type, nrow = 1) +  # scales fixed by default
  scale_x_continuous(breaks = seq(min(df_results$synergy), max(df_results$synergy), by = 0.1)) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = seq(y_min, y_max, by = 200)) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Synergy Parameter",
       y = "Total Sample Size",
       color = expression(paste(rho["AB,A"], " = ", rho["AB,B"]))) +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "gray90", color = "black"),
        legend.position = "bottom",
        legend.title = element_text(hjust = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))


# =============== Second Plot: Only FWER, with A/B/AB columns ===============
# 1) Pivot allocations into long format
df_alloc <- df_results %>%
  filter(error_type == "FWER") %>%            # keep only FWER
  dplyr::select(error_type, synergy, rho, pA_opt, pB_opt, pAB_opt) %>%
  pivot_longer(cols = c("pA_opt","pB_opt","pAB_opt"),
               names_to = "Arm",
               values_to = "Allocation")

# Rename arms
df_alloc$Arm <- factor(df_alloc$Arm,
                       levels = c("pA_opt","pB_opt","pAB_opt"),
                       labels = c("A","B","AB"))
levels(df_alloc$Arm)[levels(df_alloc$Arm) == "AB"] <- "A+B"

# 2) Plot synergy vs. Allocation, facet by Arm (A, B, AB in columns)
p_alloc <- ggplot(df_alloc, aes(x = synergy, y = Allocation,
                                color = rho, group = interaction(rho, Arm))) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = seq(min(df_alloc$synergy), max(df_alloc$synergy), by = 0.1)) +
  facet_wrap(~ Arm, nrow = 1, scales = "fixed") +  # single row, 3 columns
  scale_color_manual(values = custom_colors, drop = FALSE) +
  labs(
    x = "Synergy Parameter",
    y = "Allocation Ratio",
    color = expression(paste(rho["AB,A"], " = ", rho["AB,B"]))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "gray90", color = "black"),
    legend.position  = "bottom", legend.title = element_text(hjust = 2),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ); print(p_alloc)


