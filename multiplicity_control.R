###############################################
## 0) Libraries
###############################################
library(MASS)       # for mvrnorm
library(mvtnorm)    # for pmvnorm
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # for combining plots

###############################################
## 1) Dunnett Critical Value 
###############################################
findDunnettCritical2 <- function(rho12, alpha=0.05) {
  if (is.na(rho12) || abs(rho12) > 1) return(NA_real_)
  
  f <- function(cval) {
    if (!is.finite(cval)) return(NA_real_)
    Sigma <- matrix(c(1, rho12,
                      rho12, 1), nrow=2)
    ev <- eigen(Sigma, symmetric=TRUE)$values
    if (any(ev < -1e-12)) return(NA_real_)
    pmass <- mvtnorm::pmvnorm(
      lower=c(-cval, -cval),
      upper=c( cval,  cval),
      sigma=Sigma
    )
    if (!is.finite(pmass)) return(NA_real_)
    pmass - (1 - alpha)  # want pmass=1-alpha => difference=0
  }
  
  out <- tryCatch(
    uniroot(f, interval=c(0,5)),
    error=function(e) NULL
  )
  if (is.null(out) || !is.finite(out$root)) return(NA_real_)
  
  out$root
}

###############################################
## 2) Single-scenario simulation:
##    T1=(AB-A), T2=(B-A) with 4 methods
###############################################
simulate_arm_means <- function(nA, nB, nAB,
                               rhoAB_A, rhoAB_B, rhoA_B,
                               B=5000, alpha=0.05) {
  vA  <- 1/nA
  vB  <- 1/nB
  vAB <- 1/nAB
  
  cA_B  <- rhoA_B  / sqrt(nA*nB)
  cAB_A <- rhoAB_A / sqrt(nAB*nA)
  cAB_B <- rhoAB_B / sqrt(nAB*nB)
  
  # Cov matrix of (meanA, meanB, meanAB)
  Sigma <- matrix(c(
    vA,    cA_B,   cAB_A,
    cA_B,  vB,     cAB_B,
    cAB_A, cAB_B,  vAB
  ), nrow=3, byrow=TRUE)
  
  # Check PSD
  ev <- eigen(Sigma, symmetric=TRUE)$values
  if (any(ev < -1e-12)) {
    # invalid => return all false
    return(data.frame(
      T1=rep(NA,B), T2=rep(NA,B),
      rej1_noAdj=FALSE, rej2_noAdj=FALSE,
      rej1_bon=FALSE, rej2_bon=FALSE,
      rej1_holm=FALSE,rej2_holm=FALSE,
      rej1_dun=FALSE, rej2_dun=FALSE
    ))
  }
  
  means_mat <- MASS::mvrnorm(B, mu=c(0,0,0), Sigma=Sigma)
  meanA  <- means_mat[,1]
  meanB  <- means_mat[,2]
  meanAB <- means_mat[,3]
  
  varAB_A <- vAB + vA - 2*cAB_A
  varB_A  <- vB  + vA - 2*cA_B
  if (varAB_A<=0 || varB_A<=0) {
    # degenerate => all false
    return(data.frame(
      T1=rep(NA,B), T2=rep(NA,B),
      rej1_noAdj=FALSE, rej2_noAdj=FALSE,
      rej1_bon=FALSE, rej2_bon=FALSE,
      rej1_holm=FALSE,rej2_holm=FALSE,
      rej1_dun=FALSE, rej2_dun=FALSE
    ))
  }
  
  T1 <- (meanAB - meanA)/sqrt(varAB_A)
  T2 <- (meanB  - meanA)/sqrt(varB_A)
  
  # p-values
  p1 <- 2*pnorm(-abs(T1))
  p2 <- 2*pnorm(-abs(T2))
  
  # No adj
  rej1_noAdj <- (p1 < alpha)
  rej2_noAdj <- (p2 < alpha)
  
  # Bonf
  rej1_bon   <- (p1 < alpha/2)
  rej2_bon   <- (p2 < alpha/2)
  
  # Holm
  apply_holm_2 <- function(pp1, pp2, alpha=0.05) {
    if (is.na(pp1) || is.na(pp2)) return(c(FALSE,FALSE))
    if (pp1 < pp2) {
      r1 <- (pp1 < alpha/2)
      r2 <- if(r1) (pp2 < alpha) else FALSE
    } else {
      r2 <- (pp2 < alpha/2)
      r1 <- if(r2) (pp1 < alpha) else FALSE
    }
    c(r1,r2)
  }
  hm <- mapply(function(a,b) apply_holm_2(a,b,alpha), p1,p2)
  rej1_holm <- as.logical(hm[1,])
  rej2_holm <- as.logical(hm[2,])
  
  # Dunnett
  cov_ABA_BA <- (cAB_B - cAB_A - cA_B + vA)
  #corr_T1_T2 <- cov_ABA_BA / sqrt(varAB_A*varB_A)
  corr_T1_T2 <- 0.5
  dcval <- findDunnettCritical2(corr_T1_T2, alpha=alpha)
  if (is.na(dcval)) {
    rej1_dun <- rep(FALSE,B)
    rej2_dun <- rep(FALSE,B)
  } else {
    rej1_dun <- (abs(T1)>dcval)
    rej2_dun <- (abs(T2)>dcval)
  }
  
  data.frame(
    T1, T2,
    rej1_noAdj, rej2_noAdj,
    rej1_bon,   rej2_bon,
    rej1_holm,  rej2_holm,
    rej1_dun,   rej2_dun
  )
}

###############################################
## 3) Summarize FWER, FMER, MSFP
###############################################
summarize_fp <- function(df_out, alpha=0.05) {
  B <- nrow(df_out)
  if (B==0) {
    return(data.frame(
      fwer_noAdj=NA, fwer_bon=NA, fwer_holm=NA, fwer_dun=NA,
      fmer_noAdj=NA, fmer_bon=NA, fmer_holm=NA, fmer_dun=NA,
      msfp_noAdj=NA,msfp_bon=NA, msfp_holm=NA, msfp_dun=NA
    ))
  }
  
  # FWER
  fwer_noAdj <- mean(df_out$rej1_noAdj | df_out$rej2_noAdj, na.rm=TRUE)
  fwer_bon   <- mean(df_out$rej1_bon   | df_out$rej2_bon,   na.rm=TRUE)
  fwer_holm  <- mean(df_out$rej1_holm  | df_out$rej2_holm,  na.rm=TRUE)
  fwer_dun   <- mean(df_out$rej1_dun   | df_out$rej2_dun,   na.rm=TRUE)
  
  # FMER
  fmer_noAdj <- mean(df_out$rej1_noAdj & df_out$rej2_noAdj, na.rm=TRUE)
  fmer_bon   <- mean(df_out$rej1_bon   & df_out$rej2_bon,   na.rm=TRUE)
  fmer_holm  <- mean(df_out$rej1_holm  & df_out$rej2_holm,  na.rm=TRUE)
  fmer_dun   <- mean(df_out$rej1_dun   & df_out$rej2_dun,   na.rm=TRUE)
  
  # MSFP => "both arms declared superior in + direction"
  sup1_noAdj <- df_out$rej1_noAdj & (df_out$T1>0)
  sup2_noAdj <- df_out$rej2_noAdj & (df_out$T2>0)
  msfp_noAdj <- mean(sup1_noAdj & sup2_noAdj, na.rm=TRUE)
  
  sup1_bon <- df_out$rej1_bon & (df_out$T1>0)
  sup2_bon <- df_out$rej2_bon & (df_out$T2>0)
  msfp_bon <- mean(sup1_bon & sup2_bon, na.rm=TRUE)
  
  sup1_holm <- df_out$rej1_holm & (df_out$T1>0)
  sup2_holm <- df_out$rej2_holm & (df_out$T2>0)
  msfp_holm <- mean(sup1_holm & sup2_holm, na.rm=TRUE)
  
  sup1_dun <- df_out$rej1_dun & (df_out$T1>0)
  sup2_dun <- df_out$rej2_dun & (df_out$T2>0)
  msfp_dun <- mean(sup1_dun & sup2_dun, na.rm=TRUE)
  
  data.frame(
    fwer_noAdj, fwer_bon, fwer_holm, fwer_dun,
    fmer_noAdj, fmer_bon, fmer_holm, fmer_dun,
    msfp_noAdj, msfp_bon, msfp_holm, msfp_dun
  )
}

###############################################
## 4) Vary one correlation at a time
###############################################
investigateOneRho <- function(whichRho=c("AB_B","AB_A","A_B"),
                              rhoSeq=seq(0,0.999,by=0.05),
                              baseline1=0.3, baseline2=0.3,
                              nA=100, nB=100, nAB=100,
                              B=3000, alpha=0.05) {
  whichRho <- match.arg(whichRho)
  
  results_list <- lapply(rhoSeq, function(rVal) {
    if (whichRho=="AB_B") {
      rhoAB_B <- rVal
      rhoAB_A <- baseline1
      rhoA_B  <- baseline2
    } else if (whichRho=="AB_A") {
      rhoAB_A <- rVal
      rhoAB_B <- baseline1
      rhoA_B  <- baseline2
    } else { # "A_B"
      rhoA_B  <- rVal
      rhoAB_B <- baseline1
      rhoAB_A <- baseline2
    }
    
    df_out <- simulate_arm_means(
      nA,nB,nAB,
      rhoAB_A, rhoAB_B, 0.01,
      B=B, alpha=alpha
    )
    df_fp <- summarize_fp(df_out, alpha=alpha)
    df_fp$which_rho <- whichRho
    df_fp$rho_value <- rVal
    df_fp
  })
  
  dplyr::bind_rows(results_list)
}

###############################################
## 5) Master function: run all 3 correlations
###############################################
runAllRhos <- function(rhoSeq=seq(0,0.999,by=0.05),
                       nA=100, nB=100, nAB=100,
                       B=3000, alpha=0.05) {
  # AB_B first
  df1 <- investigateOneRho("AB_B", rhoSeq, baseline1=0.3, baseline2=0.3,
                           nA=nA,nB=nB,nAB=nAB,B=B,alpha=alpha)
  # AB_A second
  df2 <- investigateOneRho("AB_A", rhoSeq, baseline1=0.3, baseline2=0.3,
                           nA=nA,nB=nB,nAB=nAB,B=B,alpha=alpha)
  # # A_B third
  # df3 <- investigateOneRho("A_B",  rhoSeq, baseline1=0.3, baseline2=0.3,
  #                          nA=nA,nB=nB,nAB=nAB,B=B,alpha=alpha)
  bind_rows(df1,df2)
}

###############################################
## 6) run
###############################################
final_df <- runAllRhos(rhoSeq=seq(0,0.99,by=0.01),
                       nA=100, nB=100, nAB=100,
                       B=100000, alpha=0.05)

## Pivot to long format
long_df <- final_df %>%
  pivot_longer(
    cols=c("fwer_noAdj","fwer_bon","fwer_holm","fwer_dun",
           "fmer_noAdj","fmer_bon","fmer_holm","fmer_dun",
           "msfp_noAdj","msfp_bon","msfp_holm","msfp_dun"),
    names_to=c("metric","method"),
    names_pattern="(fwer|fmer|msfp)_(.*)",
    values_to="rate"
  )

## Reorder method so Bonf is drawn last (on top):
method_levels <- c("noAdj","holm","dun","bon")
method_labels <- c("NoAdj","Holm","Dunnett","Bonf") # 'Bonf' last
long_df$method <- factor(long_df$method, levels=method_levels, labels=method_labels)

## Reorder metric as FWER, FMER, MSFP
metric_levels <- c("fwer","fmer","msfp")
metric_labels <- c("FWER","FMER","MSFP")
long_df$metric <- factor(long_df$metric, levels=metric_levels, labels=metric_labels)

## Reorder which_rho as AB_B, AB_A, A_B, with math labels:
whichrho_levels <- c("AB_B","AB_A")
whichrho_labels <- c("AB_B"=expression(rho[AB*","*B]),
                     "AB_A"=expression(rho[AB*","*A]))
long_df$which_rho <- factor(long_df$which_rho, levels=whichrho_levels)


###############################################
## Final plot code 
###############################################

## 1) Define baseline rates under independent tests:
fwer_indep <- 0.0975
fmer_indep <- 0.0025
msfp_indep <- 0.000625  

## 2) Subset the data for each metric
library(dplyr)

df_fwer <- filter(long_df, metric=="FWER")
df_fmer <- filter(long_df, metric=="FMER")
df_msfp <- filter(long_df, metric=="MSFP")

# Create a shifted rho_value column only for the Holm method
df_fwer <- df_fwer %>%
  mutate(rate_shifted = ifelse(method == "Holm", rate + 0.001, rate))

## 3) Make each plot, adding geom_hline for the baseline
plot_fwer <- ggplot(df_fwer, aes(x = rho_value, y = rate_shifted, color=method)) +
  geom_line(size=1.2) +
  ## dashed line at the independent-trials   FWER
  geom_hline(yintercept = fwer_indep, linetype="dashed", color="black", size=0.9) +
  ggh4x::facet_grid2(cols = vars(which_rho), 
                     scales = "free_y", independent = "y") +  
  scale_y_continuous(limits=c(0.025,0.1)) +
  labs(x=NULL, y="FWER") +
  theme_bw(base_size=14) +
  theme(legend.position="right",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x)); plot_fwer

plot_fmer <- ggplot(df_fmer, aes(x = rho_value, y = rate, color=method)) +
  geom_line(size=1.2) +
  ## dashed line at the independent-trials FMER
  geom_hline(yintercept = fmer_indep, linetype="dashed", color="black", size=0.9) +
  ggh4x::facet_grid2(cols = vars(which_rho), 
                     scales = "free_y", independent = "y") +  
  #scale_y_continuous(limits=c(0,0.05)) +
  labs(x=NULL, y="FMER") +
  theme_bw(base_size=14) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))

plot_msfp <- ggplot(df_msfp, aes(x = rho_value, y = rate, color=method)) +
  geom_line(size=1.2) +
  ## dashed line at the independent-trials MSFP
  geom_hline(yintercept = msfp_indep, linetype="dashed", color="black", size=0.9) +
  ggh4x::facet_grid2(.~which_rho, 
                     scales = "free_y", independent = "y") + 
  #scale_y_continuous(limits=c(0,0.025)) +
  labs(x="Correlation Value", y="MSFP") +
  theme_bw(base_size=14) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))

## 4) Combine side-by-side with patchwork (or cowplot, etc.)
library(patchwork)

final_plot <- (plot_fwer | plot_fmer | plot_msfp) +
  plot_layout(nrow=3, guides="collect")

print(final_plot)
















