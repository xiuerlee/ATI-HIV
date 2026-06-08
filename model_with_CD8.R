library(deSolve)
library(tidyverse)

ati_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    T  <- max(T, 1e-12)
    L  <- max(L, 1e-12)
    I  <- max(I, 1e-12)
    V  <- max(V, 1e-12)
    CF <- max(CF, 1e-12)
    CE <- max(CE, 1e-12)
    
    K_C   <- max(K_C, 1e-12)
    K_chi <- max(K_chi, 1e-12)
    
    CD8_activation <- I / (K_C + I)
    
    CD8_exhaustion <- chi_C_G * I / (K_chi + I)
    
    CD8_killing_rate <- k_C * CF
    
    p_eff <- p_G / (1 + xi_C * CF)
    
    dT <- d_T * T0 -
      d_T * T -
      beta * V * T
    
    dL <- f_L * beta * V * T -
      d_L * L -
      a_L_G * L
    
    dI <- (1 - f_L) * beta * V * T +
      a_L_G * L -
      delta0 * I -
      CD8_killing_rate * I
    
    dV <- p_eff * I -
      c * V
    
    dCF <- s_C_G +
      rho_C_G * CD8_activation * CF -
      d_F_G * CF -
      CD8_exhaustion * CF
    
    dCE <- CD8_exhaustion * CF -
      d_E_G * CE
    
    list(c(dT, dL, dI, dV, dCF, dCE))
  })
}