library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(purrr)
library(mvtnorm)

setwd("C:/Users/mxi1/Box/Xi/abbvie/multiple testing")

# Read data
df <- read_excel("data/41591_2015_BFnm3954_MOESM10_ESM.xlsx", sheet = "PCT curve metrics")

# Process combination treatments: filter, separate and clean
df_combos <- df %>%
  # Keep only rows that contain exactly one "+"
  filter(str_count(Treatment, "\\+") == 1) %>%
  # Separate the Treatment column into two columns: DrugA, DrugB
  separate(Treatment, into = c("DrugA", "DrugB"), sep = "\\+") %>%
  # Trim any leading/trailing spaces and create the combined treatment name
  mutate(
    DrugA = str_trim(DrugA),
    DrugB = str_trim(DrugB),
    DrugA_B = paste0(DrugA, " + ", DrugB)
  ) %>%
  # Keep unique combinations
  distinct(DrugA, DrugB, DrugA_B)

# Function to compute correlations for a given trio of treatments
compute_corr <- function(drugA, drugB, drugAB, df) {
  needed_treatments <- c(drugA, drugB, drugAB)
  if (!all(needed_treatments %in% df$Treatment)) {
    return(list(
      corA_B = NA_real_,
      corA_AB = NA_real_,
      corB_AB = NA_real_
    ))
  }
  
  df_sub <- df %>%
    filter(Treatment %in% needed_treatments) %>%
    pivot_wider(
      id_cols = Model,
      names_from = Treatment,
      values_from = BestAvgResponse,
      values_fill = NA  # ensures columns are created even if missing values occur
    )
  
  if (!all(needed_treatments %in% names(df_sub))) {
    return(list(
      corA_B = NA_real_,
      corA_AB = NA_real_,
      corB_AB = NA_real_
    ))
  }
  
  df_sub <- df_sub %>%
    filter(
      !is.na(.data[[drugA]]) &
        !is.na(.data[[drugB]]) &
        !is.na(.data[[drugAB]])
    )
  
  if (nrow(df_sub) < 2) {
    return(list(
      corA_B = NA_real_,
      corA_AB = NA_real_,
      corB_AB = NA_real_
    ))
  }
  
  corA_B_val  <- cor(df_sub[[drugA]], df_sub[[drugB]])
  corA_AB_val <- cor(df_sub[[drugA]], df_sub[[drugAB]])
  corB_AB_val <- cor(df_sub[[drugB]], df_sub[[drugAB]])
  
  list(
    corA_B  = corA_B_val,
    corA_AB = corA_AB_val,
    corB_AB = corB_AB_val
  )
}

# Compute correlations for each combination row
df_combos_cor <- df_combos %>%
  rowwise() %>%
  mutate(
    corr = list(compute_corr(DrugA, DrugB, DrugA_B, df))
  ) %>%
  ungroup() %>%
  unnest_wider(corr)

# Clean up the combinations with missing or negative correlations
df_combos_cor_clean <- df_combos_cor %>%
  filter(
    !is.na(corA_B),
    !is.na(corA_AB),
    !is.na(corB_AB),
    corA_B >= 0,
    corA_AB >= 0,
    corB_AB >= 0
  )

# Count number of distinct Models for each treatment
df_combos_cor_clean <- df_combos_cor_clean %>%
  rowwise() %>%
  mutate(
    na   = df %>% filter(Treatment == DrugA) %>% distinct(Model) %>% nrow(),
    nb   = df %>% filter(Treatment == DrugB) %>% distinct(Model) %>% nrow(),
    nab  = df %>% filter(Treatment == DrugA_B) %>% distinct(Model) %>% nrow()
  ) %>%
  ungroup()

# Compute rho12 
df_combos_cor_clean <- df_combos_cor_clean %>%
  rowwise() %>%
  mutate(
    rho12 = {
      rA_B  <- corA_B
      rA_AB <- corA_AB
      rAB_B <- corB_AB
      nA    <- na
      nB    <- nb
      nAB   <- nab
      
      numerator <-
        (rAB_B / sqrt(nAB * nB)) -
        (0.001 / sqrt(nA * nB)) -
        (rA_AB / sqrt(nA * nAB)) +
        (1 / nA)
      
      term1 <- (1 / nAB) + (1 / nA) - 2 * (rAB_B / sqrt(nAB * nB))
      term2 <- (1 / nB) + (1 / nA) - 2 * (0.001 / sqrt(nA * nB))
      denominator <- sqrt(term1 * term2)
      
      if (!is.na(denominator) && denominator > 0) {
        numerator / denominator
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(rho12) & rho12 > 0 & rho12 < 1)

# --- Compute summary statistics for correlations ---
summary_stats <- df_combos_cor_clean %>%
  select(corA_B, corA_AB, corB_AB) %>%
  pivot_longer(cols = everything(), names_to = "Correlation", values_to = "Value") %>%
  group_by(Correlation) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE),
    max  = max(Value, na.rm = TRUE),
    min  = min(Value, na.rm = TRUE)
  ) %>%
  ungroup()

# For each row, compute the corresponding p-values based on rho12.
df_combos_cor_clean$p_FWER <- NA_real_
df_combos_cor_clean$p_FMER <- NA_real_
df_combos_cor_clean$p_MSFP <- NA_real_

##
## 1a) FWER (two-sided)
##
## We want c^* such that:
##   P( |Z1| <= c^*,  |Z2| <= c^* ) = 1 - alpha
## i.e. P( max(|Z1|,|Z2|) > c^* ) = alpha
## for (Z1, Z2) ~ bivariate normal with correlation rho.
##
f_cstar_fwer_2sided <- function(alpha, rho) {
  f <- function(x) {
    # Probability that both Z1 and Z2 lie in [-x, x]
    val <- pmvnorm(
      lower = c(-x, -x),
      upper = c( x,  x),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), nrow=2)
    )
    val - (1 - alpha)  # we want this = 0
  }
  # Solve numerically for x in [0, 10] (x>0 for a two-sided bound)
  uniroot(f, interval=c(0, 10))$root
}

##
## 1b) FMER (two-sided)
##
## We want c^* s.t.:
##   P( |Z1| > c^*,  |Z2| > c^* ) = alpha.
##
## That event is the union of 4 "corners":
##   (Z1>+c^* & Z2>+c^*)  OR  (Z1>+c^* & Z2<-c^*)  OR
##   (Z1<-c^* & Z2>+c^*)  OR  (Z1<-c^* & Z2<-c^*).
##
f_cstar_fmer_2sided <- function(alpha, rho) {
  f <- function(x) {
    # Probability that |Z1|>x and |Z2|>x:
    
    # corner1: Z1>=x, Z2>=x
    c1 <- pmvnorm(
      lower = c(x, x),
      upper = c(Inf, Inf),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), 2)
    )
    # corner2: Z1>=x, Z2<=-x
    c2 <- pmvnorm(
      lower = c(x, -Inf),
      upper = c(Inf, -x),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), 2)
    )
    # corner3: Z1<=-x, Z2>=x
    c3 <- pmvnorm(
      lower = c(-Inf, x),
      upper = c(-x, Inf),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), 2)
    )
    # corner4: Z1<=-x, Z2<=-x
    c4 <- pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(-x, -x),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), 2)
    )
    corners <- c1 + c2 + c3 + c4
    corners - alpha  # want = 0
  }
  uniroot(f, interval=c(0, 10))$root
}

##
## 1c) MSFP (two-sided)
##
f_cstar_msfp_2sided <- function(alpha, rho) {
  f <- function(x) {
    # Probability that |Z1|>x and |Z2|>x:
    
    # corner1: Z1>=x, Z2>=x
    c1 <- pmvnorm(
      lower = c(x, x),
      upper = c(Inf, Inf),
      mean  = c(0,0),
      corr  = matrix(c(1, rho, rho, 1), 2)
    )
    corners <- c1
    corners - alpha  # want = 0
  }
  uniroot(f, interval=c(0, 10))$root
}

for(i in seq_len(nrow(df_combos_cor_clean))) {
  rho <- df_combos_cor_clean$rho12[i]
  if(rho < 1 && rho > -1 && !is.na(rho)){
    cFWER <- f_cstar_fwer_2sided(alpha = 0.05, rho = rho)
    df_combos_cor_clean$p_FWER[i] <- 2 * (1 - pnorm(cFWER))
    
    cFMER <- f_cstar_fmer_2sided(alpha = 0.0025, rho = rho)
    df_combos_cor_clean$p_FMER[i] <- 2 * (1 - pnorm(cFMER))
    
    cMSFP <- f_cstar_msfp_2sided(alpha = 0.000625, rho = rho)
    df_combos_cor_clean$p_MSFP[i] <- 2 * (1 - pnorm(cMSFP))
  }
}

############################################################################
##         unadjusted 5% false–positive probabilities (no multiplicity
##         adjustment, common cut‑off c = qnorm(1‑0.05/2) = 1.96)
############################################################################
c_nom <- qnorm(1 - 0.05/2)      # 1.959964

prob_FWER_nom <- function(rho, c = c_nom) {
  # P( |Z1|>c OR |Z2|>c ), (Z1,Z2) bivariate N(0,0,1,1,rho)
  1 - pmvnorm(lower = c(-c, -c),
              upper = c( c,  c),
              mean  = c(0, 0),
              corr  = matrix(c(1, rho, rho, 1), 2))
}

prob_FMER_nom <- function(rho, c = c_nom) {
  # P( |Z1|>c  &  |Z2|>c )
  q <- pmvnorm(lower = c(c,  c), upper = c( Inf,  Inf),
               mean = c(0,0), corr = matrix(c(1, rho, rho, 1), 2))
  r <- pmvnorm(lower = c(c, -Inf), upper = c( Inf, -c ),
               mean = c(0,0), corr = matrix(c(1, rho, rho, 1), 2))
  s <- pmvnorm(lower = c(-Inf, c ), upper = c( -c,  Inf),
               mean = c(0,0), corr = matrix(c(1, rho, rho, 1), 2))
  t <- pmvnorm(lower = c(-Inf,-Inf), upper = c( -c, -c ),
               mean = c(0,0), corr = matrix(c(1, rho, rho, 1), 2))
  q + r + s + t                 # four “corners” of the plane
}

prob_MSFP_nom <- function(rho, c = c_nom) {
  # One‑sided multiplicity‑specific FP:
  # P( Z1>c  &  Z2>c )
  pmvnorm(lower = c(c, c), upper = c(Inf, Inf),
          mean = c(0,0),
          corr = matrix(c(1, rho, rho, 1), 2))
}

# add the three columns
df_combos_cor_clean <- df_combos_cor_clean %>%
  mutate(
    FP_FWER_nom = map_dbl(rho12,  ~ prob_FWER_nom(.x)),
    FP_FMER_nom = map_dbl(rho12,  ~ prob_FMER_nom(.x)),
    FP_MSFP_nom = map_dbl(rho12,  ~ prob_MSFP_nom(.x))
  )


# --- Compute summary statistics for single agents and combination treatments ---
# Single agents (mean and SD)
df_single <- df %>%
  filter(`Treatment type` == "single") %>%
  group_by(Treatment) %>%
  summarize(
    AvgBest = mean(BestAvgResponse, na.rm = TRUE),
    sd = sd(BestAvgResponse, na.rm = TRUE),
    n = n_distinct(Model),
    .groups = "drop"
  )

# Combinations (mean and SD)
df_combo <- df %>%
  filter(`Treatment type` == "combo") %>%
  group_by(Treatment) %>%
  summarize(
    AvgBest = mean(BestAvgResponse, na.rm = TRUE),
    sd = sd(BestAvgResponse, na.rm = TRUE),
    n = n_distinct(Model),
    .groups = "drop"
  )

# Join single agent and combination data into the combo dataset
df_combos_cor_clean <- df_combos_cor_clean %>%
  mutate(DrugAB = paste(DrugA, DrugB, sep = " + ")) %>%
  
  # Join DrugA single-agent info
  left_join(df_single, by = c("DrugA" = "Treatment")) %>%
  rename(AvgBestA = AvgBest, sdA = sd, nA = n) %>%
  
  # Join DrugB single-agent info
  left_join(df_single, by = c("DrugB" = "Treatment")) %>%
  rename(AvgBestB = AvgBest, sdB = sd, nB = n) %>%
  
  # Join combo info
  left_join(df_combo, by = c("DrugAB" = "Treatment")) %>%
  rename(AvgBestAB = AvgBest, sdAB = sd, nAB = n)

# --- Compute Pooled Standard Deviations and Normalized Effect Sizes ---
df_combos_cor_clean <- df_combos_cor_clean %>%
  rowwise() %>%
  mutate(
    # Compute pooled SD for the two single-agent treatments using:
    # pooled_sd_single = sqrt( ((nA - 1)*sdA^2 + (nB - 1)*sdB^2) / (nA + nB - 2) )
    pooled_sd_single = if_else((nA + nB - 2) > 0,
                               sqrt(((nA - 1) * sdA^2 + (nB - 1) * sdB^2) / (nA + nB - 2)),
                               NA_real_),
    
    # Compute pooled SD for comparing DrugA and the combination (DrugA + DrugB)
    pooled_sd_combo = if_else((nA + nAB - 2) > 0,
                              sqrt(((nA - 1) * sdA^2 + (nAB - 1) * sdAB^2) / (nA + nAB - 2)),
                              NA_real_),
    
    # Normalized effect sizes:
    # delta_B: difference between DrugA and DrugB normalized by the pooled SD of the two single agents
    delta_B = (AvgBestA - AvgBestB) / pooled_sd_single,
    # delta_AB: difference between DrugA and the combination (DrugA + DrugB) normalized by the pooled SD of DrugA and the combo
    delta_AB = (AvgBestA - AvgBestAB) / pooled_sd_combo,
    
    # s: ratio of the combination effect to the single agent difference
    s = delta_AB / delta_B
  ) %>%
  ungroup()

# filter out rows where the denominators or differences are not sensible
df_combos_cor_clean <- df_combos_cor_clean %>%
  filter(!is.na(delta_B), !is.na(delta_AB), delta_B > 0, delta_AB > 0)


alpha_fwer <- 0.05      # Family-Wise Error Rate (FWER) control
alpha_fmer <- 0.0025    # FMER control
alpha_msfp <- 0.000625  # MSFP control


################################################################################
## 2) Function: Probability that both |Z1| > c and |Z2| > c for bivariate normal
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

################################################################################
## 3) Function: Get critical value c* for two‐sided tests for each error type
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
      prob_both_above_c(c, rho) - alpha_msfp
    }
  } else {
    stop("Unknown error_type: must be FWER, FMER, or MSFP")
  }
  
  out <- uniroot(f_c, interval = c(0, 5))
  c_star <- out$root
  crit_cache[[key]] <<- c_star
  c_star
}

################################################################################
## 4) Function: Correlation of (Z1,Z2) under the global null
################################################################################
get_corr_Z1_Z2_null <- function(pA, pB, pAB, rhoA_B, rhoA_AB, rhoB_AB, sigma2 = 1, N = 100) {
  nA  <- round(N * pA)
  nB  <- round(N * pB)
  nAB <- round(N * pAB)
  if (min(nA, nB, nAB) < 2) return(NA)
  
  varA  <- sigma2 / nA
  varB  <- sigma2 / nB
  varAB <- sigma2 / nAB
  
  covAB_A <- rhoA_AB * sqrt(varA * varAB)
  covB_A  <- 0.001 * sqrt(varA * varB)
  covB_AB <- rhoB_AB * sqrt(varB * varAB)
  
  var_diffAB_A <- varAB + varA - 2 * covAB_A
  var_diffB_A  <- varB + varA - 2 * covB_A
  if (var_diffAB_A <= 0 || var_diffB_A <= 0) return(NA)
  
  cov_diff <- covB_AB - covB_A - covAB_A + varA
  rhoZ <- cov_diff / sqrt(var_diffAB_A * var_diffB_A)
  if (!is.finite(rhoZ) || abs(rhoZ) > 1) return(NA)
  rhoZ
}

################################################################################
## 5) Function: Simulate draws (Z1, Z2) from the alternative hypothesis
################################################################################
simulate_draws_Z <- function(pA, pB, pAB,
                             muA, muB, muAB,
                             rhoA_B, rhoA_AB, rhoB_AB, sigma2 = 1, N = 100,
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
  
  covAB_A <- rhoA_AB * sqrt(varA * varAB)
  covB_A  <- 0.001 * sqrt(varA * varB)
  covB_AB <- rhoB_AB * sqrt(varB * varAB)
  
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

################################################################################
## 6) Function: Simulate power for a given error_type using the chosen design
################################################################################
simulate_power_2sided <- function(N, pA, pB, pAB,
                                  deltaB, synergy,
                                  rhoA_B, rhoA_AB, rhoB_AB, error_type = "FWER",
                                  sigma2 = 1, nrep = 2000) {
  rhoZ <- get_corr_Z1_Z2_null(pA, pB, pAB, rhoA_B, rhoA_AB, rhoB_AB, sigma2, N);rhoZ
  if (is.na(rhoZ)) return(c(power_AB = NA, power_B = NA))
  
  c_star <- calc_crit_2sided(rhoZ, error_type = error_type);c_star
  pvalue <- 2 * (1 - pnorm(c_star))
  if(error_type == "MSFP"){
    pvalue <- 1 - pnorm(c_star)
  }
  
  muA  <- 0
  muB  <- deltaB
  muAB <- synergy * deltaB
  
  zmat <- simulate_draws_Z(pA, pB, pAB,
                           muA, muB, muAB,
                           rhoA_B, rhoA_AB, rhoB_AB, sigma2, N, nrep)
  
  valid <- complete.cases(zmat)
  if (sum(valid) == 0) return(c(power_AB = NA, power_B = NA))
  
  zmat <- zmat[valid, , drop = FALSE]
  if (error_type == "MSFP") {
    power_AB <- mean(zmat[,1] > c_star)          # 1‑sided
    power_B  <- mean(zmat[,2] > c_star)
  } else {
    power_AB <- mean(abs(zmat[,1]) > c_star)     # 2‑sided
    power_B  <- mean(abs(zmat[,2]) > c_star)
  }

  c(power_AB = power_AB, power_B = power_B, pvalue=pvalue)
}

################################################################################
## 7) Optimal Allocation via Continuous Optimization (Numerical Method)
################################################################################
optimal_allocation <- function(synergy, deltaB, rho_ABA = 0.3, rho_AB = 0.001) {
  objFun <- function(par) {
    x <- par[1]
    y <- par[2]
    pAB <- x
    pA  <- (1 - x) * y
    pB  <- 1 - pAB - pA
    # Enforce feasibility: each allocation must be at least 0.1
    if (pAB < 0.01 || pA < 0.01 || pB < 0.01) return(1e6)
    
    D1 <- (1 / pAB) + (1 / pA) - 2 * rho_ABA / sqrt(pAB * pA)
    D2 <- (1 / pB)  + (1 / pA) - 2 * 0.001   / sqrt(pA * pB)
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
evaluate_N <- function(N, synergy, rhoA_B, rhoA_AB, rhoB_AB,
                       error_type = "FWER",
                       deltaB,
                       nrep_opt = nrep_opt) {
  # Get the optimal allocation using the numerical method.
  alloc <- optimal_allocation(synergy, deltaB, rho_ABA = rhoA_AB, rho_AB = rhoA_B); alloc
  
  # Compute simulated power using these allocation probabilities.
  pow <- simulate_power_2sided(N, pA = alloc["pA"], pB = alloc["pB"], pAB = alloc["pAB"],
                               deltaB = deltaB, synergy = synergy,
                               rhoA_B=rhoA_B, rhoA_AB=rhoA_AB, rhoB_AB=rhoB_AB, error_type = error_type, nrep = nrep_opt); pow
  # Use the minimum of the two power values (for the two tests).
  power_val <- min(pow["power_AB"], pow["power_B"])
  list(
    power = power_val,
    allocation = alloc,
    pvalue=pow['pvalue']
  )
}

################################################################################
## 9) Binary search: Find minimal N that achieves target power
################################################################################
find_minN_for_80pct_binary <- function(N_start, synergy, rhoA_B, rhoA_AB, rhoB_AB,
                                       error_type = "FWER",
                                       deltaB,
                                       target_power = 0.8,
                                       nrep_opt,
                                       nrep_final,
                                       tol = -0.005) {
  # Step 1: Find an upper bound for N
  N_low <- N_start
  eval_low <- evaluate_N(N_low, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_opt);eval_low
  if (is.na(eval_low$power)) eval_low$power <- 0
  
  if (eval_low$power >= target_power) {
    N_high <- N_low
  } else {
    N_high <- N_low
    repeat {
      N_high <- N_high * 2
      eval_high <- evaluate_N(N_high, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_opt)
      if (!is.na(eval_high$power) && eval_high$power >= target_power) break
      if (N_high > 10000000) break
    }
  }
  
  # Step 2: Binary search for minimal N between N_low and N_high
  while ((N_high - N_low) > 1) {
    N_mid <- floor((N_low + N_high) / 2)
    eval_mid <- evaluate_N(N_mid, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_opt)
    if (!is.na(eval_mid$power) && eval_mid$power >= target_power) {
      N_high <- N_mid
    } else {
      N_low <- N_mid
    }
  }
  
  # Step 3: Final check with increased replication
  final_eval <- evaluate_N(N_high, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_final)
  while ((is.na(final_eval$power) || final_eval$power < target_power) && N_high < 5000) {
    N_high <- N_high + 1
    final_eval <- evaluate_N(N_high, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_final)
  }
  
  # Step 4: Try decrementing N if possible
  candidate <- N_high
  repeat {
    eval_candidate <- evaluate_N(candidate - 1, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_final)
    if (is.na(eval_candidate$power) || eval_candidate$power < target_power - tol)
      break
    candidate <- candidate - 1
    if (candidate <= N_start) break
  }
  
  final_evals <- evaluate_N(candidate, synergy, rhoA_B, rhoA_AB, rhoB_AB, error_type, deltaB, nrep_final)
  
  list(
    N = candidate,
    ratio = final_evals$allocation,
    power = final_evals$power,
    pvalue=final_evals$pvalue
  )
}

################################################################################
## 10) Main Loop: Loop over error types, synergy, and rho to determine N
################################################################################
results_list <- list()
idx <- 1

# Loop over error_type choices and each (synergy, rho) combination:
for (error_type in c("FWER", "FMER", "MSFP")) {
  for (i in seq_len(nrow(df_combos_cor_clean))) {
    s <- df_combos_cor_clean$s[i]; s
    rhoA_B <- df_combos_cor_clean$corA_B[i]; rhoA_B
    rhoA_AB <- df_combos_cor_clean$corA_AB[i]; rhoA_AB
    rhoB_AB <- df_combos_cor_clean$corB_AB[i]; rhoB_AB
    deltaB <- df_combos_cor_clean$delta_B[i]; deltaB

    #s <- 1.1
    #deltaB <- 0.4
        
    out <- find_minN_for_80pct_binary(
      N_start = 20,
      synergy = s,
      rhoA_B = rhoA_B,
      rhoA_AB = rhoA_AB,
      rhoB_AB = rhoB_AB,
      error_type = error_type,
      deltaB = deltaB,
      target_power = 0.8,
      nrep_opt = 10000,
      nrep_final = 10000,
      tol = -0.005
    ); out
    
    results_list[[idx]] <- data.frame(
      A = df_combos_cor_clean$DrugA[i],
      B = df_combos_cor_clean$DrugB[i],
      AB = df_combos_cor_clean$DrugA_B[i],
      error_type = error_type,
      delta = deltaB,
      synergy    = s,
      rhoA_B = rhoA_B,
      rhoA_AB = rhoA_AB,
      rhoB_AB = rhoB_AB,
      N          = out$N,
      pA_opt     = out$ratio["pA"],
      pB_opt     = out$ratio["pB"],
      pAB_opt    = out$ratio["pAB"],
      power      = out$power,
      pvalue=out$pvalue
    )
    idx <- idx + 1
  }
}

df_results <- do.call(rbind, results_list)

# table 2
write.csv(df_results, 'result/pdx_result.csv', row.names = F)

# table 1
write.csv(df_combos_cor_clean, 'result/pdx_fp.csv', row.names = F)

