library(mvtnorm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

###############################################
## 1) Baseline error rates under 2 independent tests
##    for alpha=0.05 (two-sided)
###############################################
fwer_indep <- 0.05  # 1 - (1-0.05)^2
fmer_indep <- 0.0025  # (0.05)^2
msfp_indep <- 0.000625 # (0.025)^2 if 2-sided "superior"

###############################################
## 2) Bivariate Probability for each metric
###############################################
clamp_rho <- function(r) {
  if (is.na(r) || r< -1) return(-1)
  if (r> 1) return(1)
  r
}

Sigma_rho <- function(rho12) {
  if (is.na(rho12) || abs(rho12)>1) return(NULL)
  matrix(c(1, rho12, rho12, 1), nrow=2)
}

fwer_prob <- function(cval, rho12) {
  rho12 <- clamp_rho(rho12)
  Sig <- Sigma_rho(rho12)
  if (is.null(Sig)) return(NA_real_)
  
  pmass <- tryCatch(
    pmvnorm(lower=c(-cval,-cval),
            upper=c(cval,cval),
            sigma=Sig),
    error=function(e) NA_real_
  )
  pmass[1]
}

fmer_prob <- function(cval, rho12) {
  rho12 <- clamp_rho(rho12)
  Sig <- Sigma_rho(rho12)
  if (is.null(Sig)) return(NA_real_)
  
  cornerProb <- function(lwr,upr) {
    tryCatch(pmvnorm(lower=lwr, upper=upr, sigma=Sig),
             error=function(e) NA_real_)
  }
  # corners: T1>c,T2>c, T1>c,T2< -c, T1< -c,T2>c, T1< -c,T2< -c
  p1 <- cornerProb(c(cval,cval), c(Inf,Inf))
  p2 <- cornerProb(c(cval,-Inf), c(Inf,-cval))
  p3 <- cornerProb(c(-Inf,cval), c(-cval,Inf))
  p4 <- cornerProb(c(-Inf,-Inf), c(-cval,-cval))
  sum(c(p1,p2,p3,p4), na.rm=FALSE)
}

msfp_prob <- function(cval, rho12) {
  rho12 <- clamp_rho(rho12)
  Sig <- Sigma_rho(rho12)
  if (is.null(Sig)) return(NA_real_)
  
  pmass <- tryCatch(
    pmvnorm(lower=c(cval,cval), upper=c(Inf,Inf), sigma=Sig),
    error=function(e) NA_real_
  )
  pmass[1]
}

###############################################
## 3) Solve for c* to match the "independent" baseline
##    e.g. FWER =>  P(|T1|<=c, |T2|<=c)= 1 - fwer_indep
###############################################
find_cstar_fwer <- function(rho12, fwer_target=0.05) {
  target <- 1 - fwer_target
  f <- function(cval) {
    val <- fwer_prob(cval, rho12)
    if (is.na(val)) return(1)
    val - target
  }
  out <- tryCatch(uniroot(f, c(0,5)), error=function(e) NULL)
  if (is.null(out)) return(NA_real_)
  out$root
}

find_cstar_fmer <- function(rho12, fmer_target=0.0025) {
  f <- function(cval) {
    val <- fmer_prob(cval, rho12)
    if (is.na(val)) return(0)
    val - fmer_target
  }
  out <- tryCatch(uniroot(f, c(0,5)), error=function(e) NULL)
  if (is.null(out)) return(NA_real_)
  out$root
}

find_cstar_msfp <- function(rho12, msfp_target=0.000625) {
  f <- function(cval) {
    val <- msfp_prob(cval, rho12)
    if (is.na(val)) return(0)
    val - msfp_target
  }
  out <- tryCatch(uniroot(f, c(0,5)), error=function(e) NULL)
  if (is.null(out)) return(NA_real_)
  out$root
}

###############################################
## 4) c* -> single-test p-value threshold
##    For FWER & FMER => 2*(1 - pnorm(c*))
##    For MSFP => 1 - pnorm(c*)
###############################################
pval_cutoff_fwer <- function(cstar) {
  if (is.na(cstar)) return(NA_real_)
  2*(1 - pnorm(cstar))
}
pval_cutoff_fmer <- function(cstar) {
  if (is.na(cstar)) return(NA_real_)
  2*(1 - pnorm(cstar))
}
pval_cutoff_msfp <- function(cstar) {
  if (is.na(cstar)) return(NA_real_)
  1 - pnorm(cstar)
}

###############################################
## 5) correlationOfT1T2
##    We'll assume nA=nB=nAB=1 for demonstration
###############################################
correlationOfT1T2 <- function(rhoAB_A, rhoAB_B, rhoA_B,
                              nA=1, nB=1, nAB=1) {
  varT1 <- (1/nAB + 1/nA - 2*rhoAB_A)
  varT2 <- (1/nB  + 1/nA - 2*rhoA_B)
  covT1T2 <- (rhoAB_B - rhoAB_A - rhoA_B + 1/nA)
  denom <- sqrt(varT1*varT2)
  if (denom<=0) return(NA_real_)
  covT1T2 / denom
}

###############################################
## 6) investigateOneRho:
##    vary one correlation param (0.01..0.99),
##    fix others at 0.3,
##    compute cstar for FWER=0.05, FMER=0.0025, MSFP=0.000625,
##    => p-value threshold
###############################################
investigateOneRho <- function(whichRho=c("AB_B","AB_A","A_B"),
                              steps=seq(0.01,0.99,by=0.01)) {
  whichRho <- match.arg(whichRho)
  
  df_list <- lapply(steps, function(rVal) {
    # fix the other 2 at 0.3
    if (whichRho=="AB_B") {
      rhoAB_B <- rVal; rhoAB_A <- 0.3; rhoA_B <- 0.01
    } else if (whichRho=="AB_A") {
      rhoAB_A <- rVal; rhoAB_B <- 0.3; rhoA_B <- 0.01
    } else {
      rhoA_B <- rVal; rhoAB_B <- 0.3; rhoAB_A <- 0.01
    }
    
    rho12 <- correlationOfT1T2(rhoAB_A, rhoAB_B, rhoA_B)
    
    cF <- find_cstar_fwer(rho12, fwer_target=0.05)
    cM <- find_cstar_fmer(rho12, fmer_target=0.0025)
    cS <- find_cstar_msfp(rho12, msfp_target=0.000625)
    
    pF <- pval_cutoff_fwer(cF)
    pM <- pval_cutoff_fmer(cM)
    pS <- pval_cutoff_msfp(cS)
    
    data.frame(
      which_rho=whichRho,
      rho_value=rVal,
      pth_fwer=pF,
      pth_fmer=pM,
      pth_msfp=pS
    )
  })
  
  do.call(rbind, df_list)
}

###############################################
## 7) runAllRhos => AB_B, AB_A, A_B
###############################################
runAllRhos <- function(steps=seq(0.01,0.99,by=0.01)) {
  df1 <- investigateOneRho("AB_B", steps=steps)
  df2 <- investigateOneRho("AB_A", steps=steps)
  #df3 <- investigateOneRho("A_B",  steps=steps)
  dplyr::bind_rows(df1, df2)
}

###############################################
## 8) usage & pivot
###############################################
final_df <- runAllRhos(steps=seq(0.01,0.99,by=0.01))

plot_df <- final_df %>%
  tidyr::pivot_longer(
    cols=c("pth_fwer","pth_fmer","pth_msfp"),
    names_to="metric",
    values_to="p_threshold"
  )

plot_df$metric <- factor(plot_df$metric,
                         levels=c("pth_fwer","pth_fmer","pth_msfp"),
                         labels=c("FWER","FMER","MSFP")
)
plot_df$which_rho <- factor(plot_df$which_rho,
                            levels=c("AB_B","AB_A"))

###############################################
## 9) Plot in 3x3 with manual y-limits
###############################################
plot_fwer <- ggplot(subset(plot_df, metric=="FWER"),
                    aes(x=rho_value, y=p_threshold)) +
  geom_line(size=1.2, color="blue") +
  facet_wrap(~which_rho, nrow=1) +
  scale_y_continuous(limits=c(0.01,0.05)) +
  labs(x=NULL, y="p-threshold\n (FWER)") +
  theme_bw(base_size=14) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))

plot_fmer <- ggplot(subset(plot_df, metric=="FMER"),
                    aes(x=rho_value, y=p_threshold)) +
  geom_line(size=1.2, color="blue") +
  facet_wrap(~which_rho, nrow=1) +
  scale_y_continuous(limits=c(0,0.05)) +
  labs(x=NULL, y="p-threshold\n (FMER)") +
  theme_bw(base_size=14) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))

plot_msfp <- ggplot(subset(plot_df, metric=="MSFP"),
                    aes(x=rho_value, y=p_threshold)) +
  geom_line(size=1.2, color="blue") +
  facet_wrap(~which_rho, nrow=1) +
  scale_y_continuous(limits=c(0,0.02)) +
  labs(x="Correlation Value", y="p-threshold\n (MSFP)") +
  theme_bw(base_size=14) +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0.05, 0.95, by = 0.15),
                     labels = function(x) sprintf("%.2f", x))

final_plot <- (plot_fwer | plot_fmer | plot_msfp) +
  plot_layout(nrow=3)

print(final_plot)
