library(deSolve)
library(tidyverse)

ati_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    T <- max(T, 1e-12)
    L <- max(L, 1e-12)
    I <- max(I, 1e-12)
    V <- max(V, 1e-12)
    
    dT <- d_T * T0 - d_T * T - beta * V * T
    
    dL <- f_L * beta * V * T -
      d_L * L -
      a_L_G * L
    
    dI <- (1 - f_L) * beta * V * T +
      a_L_G * L -
      delta0 * I -
      k_CD8 * CD8_count_G * I
    
    dV <- p_G * I -
      c * V
    
    list(c(dT, dL, dI, dV))
  })
}