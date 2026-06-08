# ============================================================
# 02_run_bayesian_full_period_model.R
# Fully Bayesian hierarchical ODE model
#
# Model:
#   Full-period infection-to-ATI model
#   CD8 cytolytic killing + CD8 non-cytolytic suppression
#   Both CD8 effects share one half-saturation constant K_C
#
# Data:
#   D:/R-4.6.0/R-ATI-project/Rdata/raw_data_from_infection.csv
#
# Output:
#   D:/R-4.6.0/R-ATI-project/bayesian_full_period_outputs
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
})

set.seed(20260607)

# ============================================================
# 0. Project directory and user settings
# ============================================================

project_dir <- "D:/R-4.6.0/R-ATI-project"
setwd(project_dir)

data_path <- "Rdata/raw_data_from_infection.csv"
out_dir <- "bayesian_full_period_outputs"

dir.create(out_dir, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Cannot find data file: ",
    data_path,
    "\nCurrent working directory is: ",
    getwd()
  )
}

# TRUE: quick test for compilation and short sampling
# FALSE: formal sampling
quick_test <- TRUE

# pVL in your CSV is original viral load, not log10 viral load.
pVL_data_is_log10 <- FALSE

# DNA in your CSV is already log10 DNA copies per 10^6 CD4 T cells.
DNA_obs_is_log10 <- TRUE

# pVL values below LOD are treated as left-censored.
use_pVL_left_censoring <- TRUE

tiny <- 1e-12

cat("[INFO] Working directory:\n")
cat("  ", getwd(), "\n")
cat("[INFO] Data path:\n")
cat("  ", data_path, "\n")

# ============================================================
# 1. Read data
# ============================================================

raw0 <- readr::read_csv(data_path, show_col_types = FALSE)

required_cols <- c(
  "id",
  "group",
  "time_dpi",
  "time_ati",
  "time_art",
  "ART_start_day",
  "ART_end_day_ATI_day0",
  "ART_on",
  "phase",
  "variable",
  "value",
  "lod"
)

missing_cols <- setdiff(required_cols, names(raw0))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

raw <- raw0 %>%
  mutate(
    id = as.character(id),
    group = tolower(as.character(group)),
    variable = as.character(variable),
    time_dpi = as.numeric(time_dpi),
    ART_start_day = as.numeric(ART_start_day),
    ART_end_day_ATI_day0 = as.numeric(ART_end_day_ATI_day0),
    value = as.numeric(value),
    lod = as.numeric(lod)
  ) %>%
  filter(
    variable %in% c("pVL", "DNA", "CD4", "CD8"),
    !is.na(id),
    !is.na(group),
    !is.na(time_dpi),
    !is.na(value),
    value > 0
  )

id_info <- raw %>%
  group_by(id) %>%
  summarise(
    group = first(group),
    ART_start_day = first(ART_start_day),
    ART_end_day_ATI_day0 = first(ART_end_day_ATI_day0),
    .groups = "drop"
  ) %>%
  arrange(group, id)

ids <- id_info$id
N_id <- length(ids)

cat("[INFO] Number of individuals:", N_id, "\n")
cat("[INFO] Observation counts:\n")
print(raw %>% count(variable))

# ============================================================
# 2. Build day-0 initial-state priors
# ============================================================

get_day0_value <- function(var_name) {
  raw %>%
    filter(variable == var_name, time_dpi == 0) %>%
    group_by(id) %>%
    summarise(value0 = first(value), .groups = "drop") %>%
    rename(!!paste0(var_name, "0") := value0)
}

get_day0_lod <- function(var_name) {
  raw %>%
    filter(variable == var_name, time_dpi == 0) %>%
    group_by(id) %>%
    summarise(lod0 = first(lod), .groups = "drop") %>%
    rename(!!paste0(var_name, "_lod0") := lod0)
}

day0_pVL <- get_day0_value("pVL")
day0_pVL_lod <- get_day0_lod("pVL")
day0_CD4 <- get_day0_value("CD4")
day0_CD8 <- get_day0_value("CD8")

init_info <- id_info %>%
  left_join(day0_pVL, by = "id") %>%
  left_join(day0_pVL_lod, by = "id") %>%
  left_join(day0_CD4, by = "id") %>%
  left_join(day0_CD8, by = "id")

median_CD4_0 <- median(init_info$CD40, na.rm = TRUE)
median_CD8_0 <- median(init_info$CD80, na.rm = TRUE)

if (is.na(median_CD4_0)) {
  stop("No day-0 CD4 values found.")
}

if (is.na(median_CD8_0)) {
  stop("No day-0 CD8 values found.")
}

init_info <- init_info %>%
  mutate(
    CD40_prior = if_else(is.na(CD40), median_CD4_0, CD40),
    CD80_prior = if_else(is.na(CD80), median_CD8_0, CD80),
    
    log_T0_prior = log(CD40_prior),
    log_C0_prior = log(CD80_prior),
    
    log_T0_prior_sd = if_else(is.na(CD40), 0.75, 0.10),
    log_C0_prior_sd = if_else(is.na(CD80), 0.75, 0.10),
    
    pVL_lod0_clean = if_else(is.na(pVL_lod0), 12.3, pVL_lod0)
  )

V0_upper <- max(init_info$pVL_lod0_clean, na.rm = TRUE)

cat("[INFO] Initial-state prior table:\n")
print(init_info)

# ============================================================
# 3. Observation model data
# ============================================================

var_levels <- c("pVL", "DNA", "CD4", "CD8")

obs_data <- raw %>%
  mutate(
    var_idx = match(variable, var_levels),
    
    is_censored = case_when(
      use_pVL_left_censoring &
        variable == "pVL" &
        !is.na(lod) &
        value < lod ~ 1L,
      TRUE ~ 0L
    ),
    
    obs_y = case_when(
      variable == "pVL" & pVL_data_is_log10 ~ value,
      variable == "pVL" & !pVL_data_is_log10 ~ log10(pmax(value, tiny)),
      
      variable == "DNA" & DNA_obs_is_log10 ~ value,
      variable == "DNA" & !DNA_obs_is_log10 ~ log10(pmax(value, tiny)),
      
      variable == "CD4" ~ log10(pmax(value, tiny)),
      variable == "CD8" ~ log10(pmax(value, tiny)),
      
      TRUE ~ NA_real_
    ),
    
    lod_y = case_when(
      variable == "pVL" & !is.na(lod) & pVL_data_is_log10 ~ lod,
      variable == "pVL" & !is.na(lod) & !pVL_data_is_log10 ~ log10(pmax(lod, tiny)),
      TRUE ~ 0
    )
  ) %>%
  filter(!is.na(obs_y))

# Day-0 pVL, CD4, CD8 are already used as initial-state priors.
# Remove them from likelihood to avoid double-counting.
obs_fit <- obs_data %>%
  filter(!(time_dpi == 0 & variable %in% c("pVL", "CD4", "CD8")))

# Global time grid:
# positive observation times + ART/ATI switch times.
time_grid <- sort(unique(c(
  obs_fit$time_dpi[obs_fit$time_dpi > 0],
  id_info$ART_start_day[id_info$ART_start_day > 0],
  id_info$ART_end_day_ATI_day0[id_info$ART_end_day_ATI_day0 > 0]
)))

N_grid <- length(time_grid)

if (N_grid < 1) {
  stop("No positive time points found.")
}

obs_fit <- obs_fit %>%
  mutate(
    id_idx = match(id, ids),
    time_idx = if_else(time_dpi == 0, 0L, match(time_dpi, time_grid))
  )

if (any(is.na(obs_fit$time_idx))) {
  stop("Some observation times were not found in time_grid.")
}

N_obs <- nrow(obs_fit)

cat("[INFO] N_grid:", N_grid, "\n")
cat("[INFO] N_obs used in likelihood:", N_obs, "\n")

# ============================================================
# 4. Prior centers
# ============================================================

fixed_pars <- list(
  d_T = 0.01,
  f_L = 0.005,
  d_L = 5e-4,
  delta0 = 1.4,
  c = 23,
  d_C = 0.02,
  eps_ART = 0.99
)

median_C <- median(raw$value[raw$variable == "CD8"], na.rm = TRUE)

if (!is.finite(median_C) || median_C <= 0) {
  median_C <- 1e6
}

p_prior_center <- 8000

pVL_proxy <- raw %>%
  filter(variable == "pVL") %>%
  mutate(
    lod_clean = if_else(is.na(lod), 1, lod)
  )

if (pVL_data_is_log10) {
  pVL_proxy <- pVL_proxy %>%
    mutate(
      V_proxy = 10^value
    )
} else {
  pVL_proxy <- pVL_proxy %>%
    mutate(
      V_proxy = pmax(value, lod_clean)
    )
}

pVL_proxy <- pVL_proxy %>%
  mutate(
    I_proxy = fixed_pars$c * V_proxy / p_prior_center
  )

median_I_proxy <- median(pVL_proxy$I_proxy, na.rm = TRUE)

if (!is.finite(median_I_proxy) || median_I_proxy <= 0) {
  median_I_proxy <- 1
}

r_T_prior_center <- 0.003

lambda_T_prior <- pmax(
  (fixed_pars$d_T - r_T_prior_center) * init_info$CD40_prior,
  1e-6
)

lambda_C_prior <- pmax(
  fixed_pars$d_C * init_info$CD80_prior,
  1e-6
)

L0_prior <- rep(1, N_id)
I0_prior <- rep(1, N_id)

# ============================================================
# 5. Stan data list
# ============================================================

stan_data <- list(
  N_id = N_id,
  N_grid = N_grid,
  N_obs = N_obs,
  
  time_grid = as.array(time_grid),
  
  obs_id = as.array(obs_fit$id_idx),
  obs_time_idx = as.array(obs_fit$time_idx),
  var_idx = as.array(obs_fit$var_idx),
  obs_y = as.array(obs_fit$obs_y),
  lod_y = as.array(obs_fit$lod_y),
  is_censored = as.array(obs_fit$is_censored),
  
  tau_ART = as.array(id_info$ART_start_day),
  tau_ATI = as.array(id_info$ART_end_day_ATI_day0),
  
  log_T0_prior = as.array(init_info$log_T0_prior),
  log_T0_prior_sd = as.array(init_info$log_T0_prior_sd),
  log_C0_prior = as.array(init_info$log_C0_prior),
  log_C0_prior_sd = as.array(init_info$log_C0_prior_sd),
  
  log_V0_upper = log(V0_upper),
  
  prior_mu_log_lambda_T = log(median(lambda_T_prior, na.rm = TRUE)),
  prior_mu_log_lambda_C = log(median(lambda_C_prior, na.rm = TRUE)),
  prior_mu_log_L0 = log(median(L0_prior, na.rm = TRUE)),
  prior_mu_log_I0 = log(median(I0_prior, na.rm = TRUE)),
  
  prior_mu_log_K_C = log(median_C),
  prior_mu_log_K_I = log(median_I_proxy),
  
  d_T = fixed_pars$d_T,
  f_L = fixed_pars$f_L,
  d_L = fixed_pars$d_L,
  delta0 = fixed_pars$delta0,
  c = fixed_pars$c,
  d_C = fixed_pars$d_C,
  eps_ART = fixed_pars$eps_ART,
  
  rel_tol = 1e-6,
  abs_tol = 1e-6,
  max_num_steps = 100000
)

# ============================================================
# 6. Write Stan model file
# ============================================================

stan_code <- '
functions {
  real log10_stan(real x) {
    return log(fmax(x, 1e-12)) / log(10.0);
  }

  vector hiv_ode(
    real t,
    vector y,
    real lambda_T,
    real lambda_C,
    real r_T,
    real beta,
    real p,
    real K_C,
    real rho_C,
    real K_I,
    real a_L,
    real k_C,
    real eta_C,
    real tau_ART,
    real tau_ATI,
    real d_T,
    real f_L,
    real d_L,
    real delta0,
    real c,
    real d_C,
    real eps_ART
  ) {
    vector[5] dydt;

    real T = fmax(y[1], 1e-12);
    real L = fmax(y[2], 1e-12);
    real I = fmax(y[3], 1e-12);
    real V = fmax(y[4], 1e-12);
    real C = fmax(y[5], 1e-12);

    real beta_eff;
    real H_C;
    real H_I;
    real p_eff;

    if (t >= tau_ART && t < tau_ATI) {
      beta_eff = beta * (1.0 - eps_ART);
    } else {
      beta_eff = beta;
    }

    H_C = C / (K_C + C);
    H_I = I / (K_I + I);

    p_eff = p * (1.0 - eta_C * H_C);
    p_eff = fmax(p_eff, 1e-12);

    dydt[1] =
      lambda_T
      + r_T * T
      - d_T * T
      - beta_eff * V * T;

    dydt[2] =
      f_L * beta_eff * V * T
      - d_L * L
      - a_L * L;

    dydt[3] =
      (1.0 - f_L) * beta_eff * V * T
      + a_L * L
      - delta0 * I
      - k_C * H_C * I;

    dydt[4] =
      p_eff * I
      - c * V;

    dydt[5] =
      lambda_C
      - d_C * C
      + rho_C * H_I * C;

    return dydt;
  }
}

data {
  int<lower=1> N_id;
  int<lower=1> N_grid;
  int<lower=1> N_obs;

  array[N_grid] real<lower=0> time_grid;

  array[N_obs] int<lower=1, upper=N_id> obs_id;
  array[N_obs] int<lower=0, upper=N_grid> obs_time_idx;
  array[N_obs] int<lower=1, upper=4> var_idx;
  array[N_obs] real obs_y;
  array[N_obs] real lod_y;
  array[N_obs] int<lower=0, upper=1> is_censored;

  array[N_id] real<lower=0> tau_ART;
  array[N_id] real<lower=0> tau_ATI;

  array[N_id] real log_T0_prior;
  array[N_id] real<lower=0> log_T0_prior_sd;
  array[N_id] real log_C0_prior;
  array[N_id] real<lower=0> log_C0_prior_sd;

  real log_V0_upper;

  real prior_mu_log_lambda_T;
  real prior_mu_log_lambda_C;
  real prior_mu_log_L0;
  real prior_mu_log_I0;

  real prior_mu_log_K_C;
  real prior_mu_log_K_I;

  real<lower=0> d_T;
  real<lower=0> f_L;
  real<lower=0> d_L;
  real<lower=0> delta0;
  real<lower=0> c;
  real<lower=0> d_C;
  real<lower=0, upper=1> eps_ART;

  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;
}

parameters {
  real<lower=log(1e-6), upper=log(0.0095)> log_r_T;
  real<lower=log(1e-12), upper=log(1e-5)> log_beta;
  real<lower=log(1), upper=log(1e7)> log_p;
  real<lower=log(1e3), upper=log(1e8)> log_K_C;
  real<lower=log(1e-8), upper=log(1)> log_rho_C;
  real<lower=log(1e-8), upper=log(1e8)> log_K_I;

  real mu_log_aL;
  real<lower=0> sigma_log_aL;
  vector[N_id] z_log_aL;

  real mu_log_kC;
  real<lower=0> sigma_log_kC;
  vector[N_id] z_log_kC;

  real mu_logit_etaC;
  real<lower=0> sigma_logit_etaC;
  vector[N_id] z_logit_etaC;

  real mu_log_lambda_T;
  real<lower=0> sigma_log_lambda_T;
  vector[N_id] z_log_lambda_T;

  real mu_log_lambda_C;
  real<lower=0> sigma_log_lambda_C;
  vector[N_id] z_log_lambda_C;

  real mu_log_L0;
  real<lower=0> sigma_log_L0;
  vector[N_id] z_log_L0;

  real mu_log_I0;
  real<lower=0> sigma_log_I0;
  vector[N_id] z_log_I0;

  vector[N_id] log_T0;
  vector[N_id] log_C0;
  vector<lower=log(1e-6), upper=log_V0_upper>[N_id] log_V0;

  vector<lower=0>[4] sigma_obs;
}

transformed parameters {
  real<lower=0> r_T;
  real<lower=0> beta;
  real<lower=0> p;
  real<lower=0> K_C;
  real<lower=0> rho_C;
  real<lower=0> K_I;

  vector<lower=0>[N_id] a_L_i;
  vector<lower=0>[N_id] k_C_i;
  vector<lower=0, upper=1>[N_id] eta_C_i;

  vector<lower=0>[N_id] lambda_T_i;
  vector<lower=0>[N_id] lambda_C_i;
  vector<lower=0>[N_id] L0_i;
  vector<lower=0>[N_id] I0_i;
  vector<lower=0>[N_id] T0_i;
  vector<lower=0>[N_id] C0_i;
  vector<lower=0>[N_id] V0_i;

  r_T = exp(log_r_T);
  beta = exp(log_beta);
  p = exp(log_p);
  K_C = exp(log_K_C);
  rho_C = exp(log_rho_C);
  K_I = exp(log_K_I);

  for (i in 1:N_id) {
    a_L_i[i] = exp(mu_log_aL + sigma_log_aL * z_log_aL[i]);
    k_C_i[i] = exp(mu_log_kC + sigma_log_kC * z_log_kC[i]);
    eta_C_i[i] = inv_logit(mu_logit_etaC + sigma_logit_etaC * z_logit_etaC[i]);

    lambda_T_i[i] = exp(mu_log_lambda_T + sigma_log_lambda_T * z_log_lambda_T[i]);
    lambda_C_i[i] = exp(mu_log_lambda_C + sigma_log_lambda_C * z_log_lambda_C[i]);

    L0_i[i] = exp(mu_log_L0 + sigma_log_L0 * z_log_L0[i]);
    I0_i[i] = exp(mu_log_I0 + sigma_log_I0 * z_log_I0[i]);

    T0_i[i] = exp(log_T0[i]);
    C0_i[i] = exp(log_C0[i]);
    V0_i[i] = exp(log_V0[i]);
  }
}

model {
  matrix[N_id * N_grid, 5] state_pred;

  // -------------------------------
  // Priors: shared parameters
  // -------------------------------

  log_r_T ~ normal(log(0.003), 0.5);
  log_beta ~ normal(log(1.58e-8), 1.0);
  log_p ~ normal(log(8000), 1.0);

  log_K_C ~ normal(prior_mu_log_K_C, 1.0);
  log_rho_C ~ normal(log(0.01), 1.0);
  log_K_I ~ normal(prior_mu_log_K_I, 1.0);

  // -------------------------------
  // Priors: hierarchical immune parameters
  // -------------------------------

  mu_log_aL ~ normal(log(1e-3), 1.0);
  sigma_log_aL ~ normal(0, 0.75);
  z_log_aL ~ std_normal();

  mu_log_kC ~ normal(log(1.0), 1.0);
  sigma_log_kC ~ normal(0, 0.75);
  z_log_kC ~ std_normal();

  mu_logit_etaC ~ normal(0, 1.5);
  sigma_logit_etaC ~ normal(0, 0.75);
  z_logit_etaC ~ std_normal();

  // -------------------------------
  // Priors: nuisance parameters
  // -------------------------------

  mu_log_lambda_T ~ normal(prior_mu_log_lambda_T, 1.5);
  sigma_log_lambda_T ~ normal(0, 0.75);
  z_log_lambda_T ~ std_normal();

  mu_log_lambda_C ~ normal(prior_mu_log_lambda_C, 1.5);
  sigma_log_lambda_C ~ normal(0, 0.75);
  z_log_lambda_C ~ std_normal();

  mu_log_L0 ~ normal(prior_mu_log_L0, 2.0);
  sigma_log_L0 ~ normal(0, 1.0);
  z_log_L0 ~ std_normal();

  mu_log_I0 ~ normal(prior_mu_log_I0, 2.0);
  sigma_log_I0 ~ normal(0, 1.0);
  z_log_I0 ~ std_normal();

  for (i in 1:N_id) {
    log_T0[i] ~ normal(log_T0_prior[i], log_T0_prior_sd[i]);
    log_C0[i] ~ normal(log_C0_prior[i], log_C0_prior_sd[i]);

    // Day-0 pVL is often below LOD.
    // Here V0 is weakly estimated below the LOD upper bound.
    log_V0[i] ~ normal(log(1.0), 1.0);
  }

  // residual SDs on log10 observation scale
  sigma_obs[1] ~ normal(0, 0.75);  // pVL
  sigma_obs[2] ~ normal(0, 0.50);  // DNA
  sigma_obs[3] ~ normal(0, 0.25);  // CD4
  sigma_obs[4] ~ normal(0, 0.25);  // CD8

  // -------------------------------
  // ODE solve for each individual
  // -------------------------------

  for (i in 1:N_id) {
    vector[5] y0;
    array[N_grid] vector[5] sol;

    y0[1] = T0_i[i];
    y0[2] = L0_i[i];
    y0[3] = I0_i[i];
    y0[4] = V0_i[i];
    y0[5] = C0_i[i];

    sol = ode_bdf_tol(
      hiv_ode,
      y0,
      0.0,
      time_grid,
      rel_tol,
      abs_tol,
      max_num_steps,
      lambda_T_i[i],
      lambda_C_i[i],
      r_T,
      beta,
      p,
      K_C,
      rho_C,
      K_I,
      a_L_i[i],
      k_C_i[i],
      eta_C_i[i],
      tau_ART[i],
      tau_ATI[i],
      d_T,
      f_L,
      d_L,
      delta0,
      c,
      d_C,
      eps_ART
    );

    for (j in 1:N_grid) {
      int row_id = (i - 1) * N_grid + j;

      state_pred[row_id, 1] = fmax(sol[j][1], 1e-12);
      state_pred[row_id, 2] = fmax(sol[j][2], 1e-12);
      state_pred[row_id, 3] = fmax(sol[j][3], 1e-12);
      state_pred[row_id, 4] = fmax(sol[j][4], 1e-12);
      state_pred[row_id, 5] = fmax(sol[j][5], 1e-12);
    }
  }

  // -------------------------------
  // Observation likelihood
  //
  // variable index:
  // 1 = pVL
  // 2 = DNA
  // 3 = CD4
  // 4 = CD8
  // -------------------------------

  for (n in 1:N_obs) {
    int i = obs_id[n];
    int tidx = obs_time_idx[n];

    real T;
    real L;
    real I;
    real V;
    real C;

    real pred;
    real denom;

    if (tidx == 0) {
      T = T0_i[i];
      L = L0_i[i];
      I = I0_i[i];
      V = V0_i[i];
      C = C0_i[i];
    } else {
      int row_id = (i - 1) * N_grid + tidx;

      T = state_pred[row_id, 1];
      L = state_pred[row_id, 2];
      I = state_pred[row_id, 3];
      V = state_pred[row_id, 4];
      C = state_pred[row_id, 5];
    }

    denom = fmax(T + L + I, 1e-12);

    if (var_idx[n] == 1) {
      pred = log10_stan(V);
    } else if (var_idx[n] == 2) {
      pred = log10_stan(1e6 * (L + I) / denom);
    } else if (var_idx[n] == 3) {
      pred = log10_stan(T + L + I);
    } else {
      pred = log10_stan(C);
    }

    if (is_censored[n] == 1) {
      target += normal_lcdf(lod_y[n] | pred, sigma_obs[var_idx[n]]);
    } else {
      target += normal_lpdf(obs_y[n] | pred, sigma_obs[var_idx[n]]);
    }
  }
}
'

stan_file <- file.path(out_dir, "bayesian_full_period_cd8_suppression.stan")
writeLines(stan_code, stan_file)

cat("[INFO] Stan model written to:\n")
cat("  ", stan_file, "\n")

# ============================================================
# 7. Compile Stan model
# ============================================================

mod <- cmdstan_model(stan_file)

# ============================================================
# 8. Initial values
# ============================================================

init_fun <- function(chain_id = 1) {
  list(
    log_r_T = log(0.003),
    log_beta = log(1.58e-8),
    log_p = log(8000),
    log_K_C = log(median_C),
    log_rho_C = log(0.01),
    log_K_I = log(median_I_proxy),
    
    mu_log_aL = log(1e-3),
    sigma_log_aL = 0.2,
    z_log_aL = rep(0, N_id),
    
    mu_log_kC = log(1),
    sigma_log_kC = 0.2,
    z_log_kC = rep(0, N_id),
    
    mu_logit_etaC = 0,
    sigma_logit_etaC = 0.2,
    z_logit_etaC = rep(0, N_id),
    
    mu_log_lambda_T = stan_data$prior_mu_log_lambda_T,
    sigma_log_lambda_T = 0.2,
    z_log_lambda_T = rep(0, N_id),
    
    mu_log_lambda_C = stan_data$prior_mu_log_lambda_C,
    sigma_log_lambda_C = 0.2,
    z_log_lambda_C = rep(0, N_id),
    
    mu_log_L0 = stan_data$prior_mu_log_L0,
    sigma_log_L0 = 0.5,
    z_log_L0 = rep(0, N_id),
    
    mu_log_I0 = stan_data$prior_mu_log_I0,
    sigma_log_I0 = 0.5,
    z_log_I0 = rep(0, N_id),
    
    log_T0 = init_info$log_T0_prior,
    log_C0 = init_info$log_C0_prior,
    log_V0 = rep(log(1), N_id),
    
    sigma_obs = c(0.5, 0.3, 0.1, 0.1)
  )
}

# ============================================================
# 9. Sampling
# ============================================================

if (quick_test) {
  chains <- 2
  parallel_chains <- 2
  iter_warmup <- 200
  iter_sampling <- 200
  adapt_delta <- 0.90
  max_treedepth <- 11
} else {
  chains <- 4
  parallel_chains <- 4
  iter_warmup <- 1000
  iter_sampling <- 1000
  adapt_delta <- 0.95
  max_treedepth <- 13
}

fit <- mod$sample(
  data = stan_data,
  seed = 20260607,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  init = init_fun,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = 50
)

# ============================================================
# 10. Save diagnostics and posterior summaries
# ============================================================

fit$save_object(file = file.path(out_dir, "fit_bayesian_full_period.RDS"))

summary_all <- fit$summary()

readr::write_csv(
  summary_all,
  file.path(out_dir, "posterior_summary_all_parameters.csv")
)

core_vars <- c(
  "r_T",
  "beta",
  "p",
  "K_C",
  "rho_C",
  "K_I",
  "mu_log_aL",
  "sigma_log_aL",
  "mu_log_kC",
  "sigma_log_kC",
  "mu_logit_etaC",
  "sigma_logit_etaC",
  "sigma_obs"
)

summary_core <- fit$summary(variables = core_vars)

readr::write_csv(
  summary_core,
  file.path(out_dir, "posterior_summary_core_population_parameters.csv")
)

summary_individual <- fit$summary(
  variables = c("a_L_i", "k_C_i", "eta_C_i")
)

readr::write_csv(
  summary_individual,
  file.path(out_dir, "posterior_summary_individual_immune_parameters.csv")
)

# ============================================================
# 11. Posterior early/late comparison
# ============================================================

draws_ind <- fit$draws(
  variables = c("a_L_i", "k_C_i", "eta_C_i"),
  format = "df"
) %>%
  as.data.frame()

early_idx <- which(id_info$group == "early")
late_idx <- which(id_info$group == "late")

extract_matrix <- function(draws_df, prefix, N_id) {
  cols <- paste0(prefix, "[", seq_len(N_id), "]")
  missing <- setdiff(cols, names(draws_df))
  
  if (length(missing) > 0) {
    stop("Missing posterior columns: ", paste(missing, collapse = ", "))
  }
  
  as.matrix(draws_df[, cols, drop = FALSE])
}

summarise_individual_param <- function(mat, parameter_name) {
  tibble(
    id = ids,
    group = id_info$group,
    parameter = parameter_name,
    median = apply(mat, 2, median),
    q025 = apply(mat, 2, quantile, probs = 0.025),
    q975 = apply(mat, 2, quantile, probs = 0.975)
  )
}

summarise_group_comparison <- function(mat, parameter_name) {
  early_mean <- rowMeans(mat[, early_idx, drop = FALSE])
  late_mean <- rowMeans(mat[, late_idx, drop = FALSE])
  ratio <- early_mean / late_mean
  diff <- early_mean - late_mean
  
  tibble(
    parameter = parameter_name,
    
    early_median = median(early_mean),
    early_q025 = quantile(early_mean, 0.025),
    early_q975 = quantile(early_mean, 0.975),
    
    late_median = median(late_mean),
    late_q025 = quantile(late_mean, 0.025),
    late_q975 = quantile(late_mean, 0.975),
    
    ratio_early_over_late_median = median(ratio),
    ratio_early_over_late_q025 = quantile(ratio, 0.025),
    ratio_early_over_late_q975 = quantile(ratio, 0.975),
    
    diff_early_minus_late_median = median(diff),
    diff_early_minus_late_q025 = quantile(diff, 0.025),
    diff_early_minus_late_q975 = quantile(diff, 0.975),
    
    posterior_prob_early_gt_late = mean(early_mean > late_mean),
    posterior_prob_early_lt_late = mean(early_mean < late_mean)
  )
}

mat_aL <- extract_matrix(draws_ind, "a_L_i", N_id)
mat_kC <- extract_matrix(draws_ind, "k_C_i", N_id)
mat_eta <- extract_matrix(draws_ind, "eta_C_i", N_id)

individual_param_summary <- bind_rows(
  summarise_individual_param(mat_aL, "a_L"),
  summarise_individual_param(mat_kC, "k_C"),
  summarise_individual_param(mat_eta, "eta_C")
)

group_comparison_summary <- bind_rows(
  summarise_group_comparison(mat_aL, "a_L"),
  summarise_group_comparison(mat_kC, "k_C"),
  summarise_group_comparison(mat_eta, "eta_C")
)

readr::write_csv(
  individual_param_summary,
  file.path(out_dir, "posterior_individual_immune_parameter_intervals.csv")
)

readr::write_csv(
  group_comparison_summary,
  file.path(out_dir, "posterior_early_late_comparison.csv")
)

# ============================================================
# 12. Basic diagnostic print
# ============================================================

cat("\n[DONE] Bayesian model finished.\n")
cat("[DONE] Outputs written to:\n")
cat("  ", out_dir, "\n\n")

cat("[IMPORTANT] Check these files first:\n")
cat("  1. posterior_summary_core_population_parameters.csv\n")
cat("  2. posterior_summary_individual_immune_parameters.csv\n")
cat("  3. posterior_individual_immune_parameter_intervals.csv\n")
cat("  4. posterior_early_late_comparison.csv\n")
cat("  5. posterior_summary_all_parameters.csv\n\n")

cat("[IMPORTANT] Check convergence columns:\n")
cat("  - rhat should be close to 1.00\n")
cat("  - ess_bulk and ess_tail should not be too small\n")
cat("  - if there are divergent transitions, increase adapt_delta or simplify model\n\n")

print(group_comparison_summary)