library(deSolve)
library(tidyverse)

ati_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    T <- max(T, 1e-12)
    L <- max(L, 1e-12)
    I <- max(I, 1e-12)
    V <- max(V, 1e-12)
    
    CD8_count_G <- max(CD8_count_G, 1e-12)
    
    CD107_percent_G <- max(CD107_percent_G, 0)
    IFNg_percent_G <- max(IFNg_percent_G, 0)
    
    CD107_frac_G <- CD107_percent_G / 100
    IFNg_frac_G <- IFNg_percent_G / 100
    
    CD107_frac_G <- min(CD107_frac_G, 1)
    IFNg_frac_G <- min(IFNg_frac_G, 1)
    
    K_CD107 <- max(K_CD107, 1e-12)
    K_IFNg <- max(K_IFNg, 1e-12)
    
    CD107_effector <- CD8_count_G * CD107_frac_G
    IFNg_effector <- CD8_count_G * IFNg_frac_G
    
    CD8_killing_rate <- k_CD8 *
      CD107_effector / (K_CD107 + CD107_effector)
    
    IFNg_suppression <- IFNg_effector / (K_IFNg + IFNg_effector)
    
    p_eff <- p_G / (1 + eta_IFNg * IFNg_suppression)
    
    dT <- d_T * T0 - d_T * T - beta * V * T
    
    dL <- f_L * beta * V * T -
      d_L * L -
      a_L_G * L
    
    dI <- (1 - f_L) * beta * V * T +
      a_L_G * L -
      delta0 * I -
      CD8_killing_rate * I
    
    dV <- p_eff * I -
      c * V
    
    list(c(dT, dL, dI, dV))
  })
}