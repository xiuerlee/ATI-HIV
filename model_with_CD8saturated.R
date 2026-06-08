library(deSolve)
library(tidyverse)

ati_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    T <- max(T, 1e-12)
    L <- max(L, 1e-12)
    I <- max(I, 1e-12)
    V <- max(V, 1e-12)
    
    CD8_count_G <- max(CD8_count_G, 1e-12)
    K_CD8 <- max(K_CD8, 1e-12)
    
    CD8_killing_rate <- k_CD8 * CD8_count_G / (K_CD8 + CD8_count_G)
    
    dT <- d_T * T0 - d_T * T - beta * V * T
    
    dL <- f_L * beta * V * T -
      d_L * L -
      a_L_G * L
    
    dI <- (1 - f_L) * beta * V * T +
      a_L_G * L -
      delta0 * I -
      CD8_killing_rate * I
    
    dV <- p_G * I -
      c * V
    
    list(c(dT, dL, dI, dV))
  })
}