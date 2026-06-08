library(deSolve)
library(tidyverse)

source("model_without_CD8.R")

data_clean <- read.csv("../Rdata/raw_data.csv")
init_table <- read.csv("../Rdata/init_table.csv")

fixed_pars <- c(
  d_T = 0.01,
  beta = 1.58e-8,
  f_L = 0.005,
  d_L = 5e-4,
  c = 23
)

test_pars <- c(
  fixed_pars,
  a_L = 1e-3,
  delta = 1.4,
  p = 6e3,
  theta_a = 0,
  theta_delta = 0,
  theta_p = 0
)

test_flags <- c(
  za = 0,
  zdelta = 0,
  zp = 0
)

simulate_one <- function(id_row, pars, times, model_flags) {
  
  G <- ifelse(id_row$group == "early", 1, 0)
  
  a_L_G <- as.numeric(pars["a_L"]) *
    exp(-as.numeric(model_flags["za"]) * as.numeric(pars["theta_a"]) * G)
  
  delta_G <- as.numeric(pars["delta"]) *
    exp(as.numeric(model_flags["zdelta"]) * as.numeric(pars["theta_delta"]) * G)
  
  p_G <- as.numeric(pars["p"]) *
    exp(-as.numeric(model_flags["zp"]) * as.numeric(pars["theta_p"]) * G)
  
  V0 <- max(as.numeric(id_row$V0), 1e-12)
  L0 <- max(as.numeric(id_row$L0), 1e-12)
  T0 <- max(as.numeric(id_row$T0), 1e-12)
  
  I0 <- as.numeric(pars["c"]) * V0 / p_G
  I0 <- max(I0, 1e-12)
  
  state <- c(
    T = T0,
    L = L0,
    I = I0,
    V = V0
  )
  
  pars_use <- c(
    pars,
    T0 = T0,
    a_L_G = a_L_G,
    delta_G = delta_G,
    p_G = p_G
  )
  
  out <- tryCatch({
    ode(
      y = state,
      times = times,
      func = ati_model,
      parms = pars_use,
      method = "lsoda"
    )
  }, error = function(e) {
    message("ODE error: ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) {
    return(NULL)
  }
  
  as.data.frame(out) %>%
    mutate(
      id = id_row$id,
      group = id_row$group,
      a_L_G = a_L_G,
      delta_G = delta_G,
      p_G = p_G,
      CD4_total = T + L + I,
      pred_pVL_log10 = log10(pmax(V, 1e-12)),
      pred_DNA_log10 = log10(pmax(1e6 * (L + I) / CD4_total, 1e-12)),
      pred_CD4_log10 = log10(pmax(CD4_total, 1e-12))
    )
}

times <- seq(0, max(data_clean$time), by = 0.1)

sim_test <- simulate_one(
  id_row = init_table[1, ],
  pars = test_pars,
  times = times,
  model_flags = test_flags
)

if (is.null(sim_test)) {
  stop("Simulation failed.")
}

sim_test %>%
  pivot_longer(
    cols = c(T, L, I, V, CD4_total),
    names_to = "state",
    values_to = "value"
  ) %>%
  ggplot(aes(time, log10(pmax(value, 1e-12)), color = state)) +
  geom_line() +
  theme_bw()
