library(deSolve)
library(tidyverse)
library(minpack.lm)

# ============================================================
# Total-CD8 saturated-killing fixed-parameter calibration script
# with xi_C non-cytolytic suppression retained
# ============================================================
# Purpose:
#   This script is NOT for model selection.
#   It tests several biologically plausible fixed-parameter settings
#   before later model selection.
#
# Current model:
#   T, L, I, V, C
#   C = total CD8 state variable
#
# CD8 effects:
#   1) Cytolytic killing of productively infected cells:
#        k_C * (C_scaled / (K_kill + C_scaled)) * I
#
#   2) Non-cytolytic suppression of viral production:
#        p_eff = p / (1 + xi_C * C_scaled)
#
#   where:
#        C_scaled = C / 1e6
#
# Important:
#   This script fits common parameters only.
#   It does NOT test early/late group differences.
#   It is meant to choose reasonable fixed parameters before model selection.
#
# Fitted common parameters in each fixed-parameter setting:
#   a_L, p, k_C, xi_C, rho_C
#
# Fixed parameters tested across settings:
#   delta0, c, d_C, K_C, K_kill
# ============================================================

# ============================================================
# 0. Global settings
# ============================================================

EPS <- 1e-12

N_RANDOM_START <- 5
RANDOM_SEED_BASE <- 1000

USE_VARIABLE_SCALING <- TRUE
FIT_VARIABLES <- c("pVL", "DNA", "CD4", "CD8")

# Keep this FALSE for the first calibration pass.
# If TRUE, I0 is recalculated as c*V0/p for each fitted p and fixed c.
# This changes the initial condition and can confound DNA fitting.
RECOMPUTE_I0_FROM_QSS <- FALSE

CD8_SCALE <- 1e6

# Current-directory workflow:
# Put this script, raw_data_with_cd8_counts.csv, and
# init_table_total_CD8_no_exhaustion.csv in the same folder.
DATA_DIR <- "D:/R-4.6.0/R-ATI-project/Rdata"

DATA_FILE <- file.path(DATA_DIR, "raw_data_with_cd8_counts.csv")
INIT_FILE <- file.path(DATA_DIR, "init_table_total_CD8_no_exhaustion.csv")

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  DATA_DIR,
  "fixed_parameter_calibration_total_CD8_saturated_killing_with_xiC",
  run_tag
)

fig_dir <- file.path(out_dir, "figures")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Fixed-parameter sets to test
# ============================================================
# Interpretation:
#   delta0 : natural loss rate of productively infected cells, day^-1
#   c      : free virion clearance rate, day^-1
#   d_C    : CD8 return-to-baseline/homeostatic rate, day^-1
#   K_C    : half-saturation for antigen-driven CD8 activation; same internal scale as I
#   K_kill : half-saturation for CD8 killing; C_scaled = C / 1e6
#
# Why these sets:
#   - pVL was previously underpredicted. Too large delta0, too strong k_C,
#     too small K_kill, or too large xi_C can all push V downward.
#   - CD8 was previously overpredicted. Too small K_C makes CD8 activation
#     nearly saturated whenever I > K_C. The old K_C = 100 may be too small
#     if I is on a cells/mL-like internal scale.
#   - K_kill > 1 weakens effective killing saturation when C_scaled is around 1.

fixed_parameter_grid <- tribble(
  ~setting_id, ~description,                 ~d_T,  ~beta,     ~f_L,   ~d_L,   ~delta0, ~c,  ~d_C, ~K_C,  ~K_kill, ~CD8_scale,
  "S01",      "original_reference",          0.01, 1.58e-8,   0.005,  5e-4,   1.40,    23,  0.75, 1e2,   1.0,     CD8_SCALE,
  "S02",      "lower_delta_reduce_CD8_act",  0.01, 1.58e-8,   0.005,  5e-4,   0.70,    23,  0.75, 1e4,   5.0,     CD8_SCALE,
  "S03",      "classic_delta_reduce_kill",   0.01, 1.58e-8,   0.005,  5e-4,   0.45,    23,  0.75, 1e4,   5.0,     CD8_SCALE,
  "S04",      "moderate_c_lower_delta",      0.01, 1.58e-8,   0.005,  5e-4,   0.70,    10,  0.75, 1e4,   5.0,     CD8_SCALE,
  "S05",      "low_c_classic_delta",         0.01, 1.58e-8,   0.005,  5e-4,   0.45,    3,   0.75, 1e4,   5.0,     CD8_SCALE,
  "S06",      "reduce_CD8_activation",       0.01, 1.58e-8,   0.005,  5e-4,   1.00,    23,  0.75, 1e5,   5.0,     CD8_SCALE,
  "S07",      "reduce_CD8_killing",          0.01, 1.58e-8,   0.005,  5e-4,   1.00,    23,  0.75, 1e4,   10.0,    CD8_SCALE,
  "S08",      "reduce_act_and_kill",         0.01, 1.58e-8,   0.005,  5e-4,   0.70,    23,  1.00, 1e5,   10.0,    CD8_SCALE,
  "S09",      "faster_CD8_homeostasis",      0.01, 1.58e-8,   0.005,  5e-4,   0.70,    23,  1.00, 1e4,   5.0,     CD8_SCALE,
  "S10",      "slower_CD8_homeostasis",      0.01, 1.58e-8,   0.005,  5e-4,   0.70,    23,  0.25, 1e4,   5.0,     CD8_SCALE,
  "S11",      "balanced_mid",                0.01, 1.58e-8,   0.005,  5e-4,   0.70,    10,  1.00, 1e5,   10.0,    CD8_SCALE,
  "S12",      "moderate_all",                0.01, 1.58e-8,   0.005,  5e-4,   1.00,    10,  0.75, 1e4,   5.0,     CD8_SCALE
)

write.csv(
  fixed_parameter_grid,
  file.path(out_dir, "fixed_parameter_test_grid.csv"),
  row.names = FALSE
)

# ============================================================
# 2. Read data
# ============================================================
# Expected raw data long format:
# id, group, time, variable, value, lod
#
# pVL: raw copies/mL
# DNA: already log10 copies per 10^6 CD4 cells
# CD4: raw cells/mL
# CD8: raw cells/mL

if (!file.exists(DATA_FILE)) {
  stop("找不到 DATA_FILE: ", DATA_FILE, "\n当前工作目录是: ", getwd())
}

if (!file.exists(INIT_FILE)) {
  stop("找不到 INIT_FILE: ", INIT_FILE, "\n当前工作目录是: ", getwd())
}

data_clean <- read.csv(DATA_FILE)
init_table <- read.csv(INIT_FILE)

# ============================================================
# 3. Helper functions
# ============================================================

is_early_group <- function(group_value) {
  tolower(as.character(group_value)) %in% c(
    "early",
    "early_art",
    "early art",
    "early-treated",
    "early treated",
    "et",
    "w4"
  )
}

safe_log10 <- function(x) {
  log10(pmax(x, EPS))
}

save_plot_png_pdf <- function(plot_obj, filename_base, width = 8, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(fig_dir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height
  )
}

row_to_fixed_pars <- function(parameter_row) {
  c(
    d_T = as.numeric(parameter_row$d_T),
    beta = as.numeric(parameter_row$beta),
    f_L = as.numeric(parameter_row$f_L),
    d_L = as.numeric(parameter_row$d_L),
    delta0 = as.numeric(parameter_row$delta0),
    c = as.numeric(parameter_row$c),
    d_C = as.numeric(parameter_row$d_C),
    K_C = as.numeric(parameter_row$K_C),
    K_kill = as.numeric(parameter_row$K_kill),
    CD8_scale = as.numeric(parameter_row$CD8_scale)
  )
}

# ============================================================
# 4. Data cleaning
# ============================================================

data_raw <- data_clean %>%
  mutate(
    id = as.character(id),
    group = as.character(group),
    time = as.numeric(time),
    variable = trimws(as.character(variable)),
    value = as.numeric(value),
    lod = as.numeric(lod)
  )

init_raw <- init_table %>%
  mutate(
    id = as.character(id),
    group = as.character(group),
    T0 = as.numeric(T0),
    L0 = as.numeric(L0),
    I0 = as.numeric(I0),
    V0 = as.numeric(V0),
    CD8_0 = if ("CD8_0" %in% names(.)) as.numeric(CD8_0) else NA_real_,
    CD8_0_source = if ("CD8_0_source" %in% names(.)) as.character(CD8_0_source) else NA_character_
  )

init_excluded <- init_raw %>%
  filter(is.na(T0) | is.na(L0) | is.na(I0) | is.na(V0))

if (nrow(init_excluded) > 0) {
  write.csv(
    init_excluded,
    file.path(out_dir, "init_table_excluded_due_to_missing_initials.csv"),
    row.names = FALSE
  )
  
  warning(
    "Some individuals have missing T0/L0/I0/V0 and will be excluded. See init_table_excluded_due_to_missing_initials.csv."
  )
}

init_raw <- init_raw %>%
  filter(!is.na(T0), !is.na(L0), !is.na(I0), !is.na(V0))

# ============================================================
# 5. Data scale checks
# ============================================================

unit_check <- data_raw %>%
  filter(variable %in% c("pVL", "DNA", "CD4", "CD8")) %>%
  group_by(variable, group) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(unit_check)

write.csv(
  unit_check,
  file.path(out_dir, "input_unit_check.csv"),
  row.names = FALSE
)

if (
  any(
    unit_check$variable == "CD4" &
    is.finite(unit_check$median_value) &
    unit_check$median_value < 10000
  )
) {
  warning(
    "CD4 median is < 10000. CD4 may still be in cells/uL rather than cells/mL. ",
    "The current script assumes CD4 is already cells/mL."
  )
}

if (
  any(
    unit_check$variable == "CD8" &
    is.finite(unit_check$median_value) &
    unit_check$median_value < 10000
  )
) {
  warning(
    "CD8 median is < 10000. CD8 may still be in cells/uL rather than cells/mL. ",
    "The current script assumes CD8 is already cells/mL."
  )
}

# ============================================================
# 6. Add initial CD8 to initial table
# ============================================================
# Priority:
#   1) CD8_0 already provided in init table
#   2) individual CD8 at time 0
#   3) individual median CD8
#   4) group median CD8
#   5) overall median CD8

cd8_data <- data_raw %>%
  filter(variable == "CD8", !is.na(value), is.finite(value))

if (nrow(cd8_data) == 0) {
  stop("No CD8 data found in raw_data_with_cd8_counts.csv.")
}

overall_cd8_median <- median(cd8_data$value, na.rm = TRUE)

cd8_by_group <- cd8_data %>%
  group_by(group) %>%
  summarise(
    CD8_count_group = median(value, na.rm = TRUE),
    .groups = "drop"
  )

cd8_by_id <- cd8_data %>%
  group_by(id, group) %>%
  summarise(
    CD8_time0 = ifelse(
      any(time == 0, na.rm = TRUE),
      median(value[time == 0], na.rm = TRUE),
      NA_real_
    ),
    CD8_median_id = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CD8_0_source_id = case_when(
      is.finite(CD8_time0) ~ "time0",
      is.finite(CD8_median_id) ~ "id_median",
      TRUE ~ NA_character_
    ),
    CD8_count_id = case_when(
      is.finite(CD8_time0) ~ CD8_time0,
      is.finite(CD8_median_id) ~ CD8_median_id,
      TRUE ~ NA_real_
    )
  ) %>%
  select(id, group, CD8_count_id, CD8_0_source_id)

init_fit <- init_raw %>%
  left_join(cd8_by_id, by = c("id", "group")) %>%
  left_join(cd8_by_group, by = "group") %>%
  mutate(
    C0 = case_when(
      is.finite(CD8_0) ~ CD8_0,
      is.finite(CD8_count_id) ~ CD8_count_id,
      is.finite(CD8_count_group) ~ CD8_count_group,
      TRUE ~ overall_cd8_median
    ),
    C0_source = case_when(
      is.finite(CD8_0) & !is.na(CD8_0_source) ~ CD8_0_source,
      is.finite(CD8_0) ~ "init_table_CD8_0",
      !is.finite(CD8_0) & is.finite(CD8_count_id) ~ CD8_0_source_id,
      !is.finite(CD8_0) & !is.finite(CD8_count_id) & is.finite(CD8_count_group) ~ "group_median",
      TRUE ~ "overall_median"
    )
  ) %>%
  filter(is.finite(C0), C0 > 0)

write.csv(
  init_fit,
  file.path(out_dir, "init_table_with_total_CD8_used.csv"),
  row.names = FALSE
)

valid_ids <- init_fit$id

# ============================================================
# 7. Fitting data preprocessing
# ============================================================

data_fit <- data_raw %>%
  filter(
    id %in% valid_ids,
    variable %in% FIT_VARIABLES,
    !is.na(time),
    !is.na(value),
    is.finite(time),
    is.finite(value)
  ) %>%
  mutate(
    value_raw = value,
    lod_raw = lod,
    
    value = case_when(
      variable == "pVL" ~ safe_log10(value_raw),
      variable == "DNA" ~ value_raw,
      variable == "CD4" ~ safe_log10(value_raw),
      variable == "CD8" ~ safe_log10(value_raw),
      TRUE ~ value_raw
    ),
    
    lod = case_when(
      variable == "pVL" &
        !is.na(lod_raw) &
        is.finite(lod_raw) ~ safe_log10(lod_raw),
      TRUE ~ NA_real_
    ),
    
    is_censored = case_when(
      variable == "pVL" &
        !is.na(lod_raw) &
        is.finite(lod_raw) &
        value_raw <= lod_raw + 1e-10 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is.na(value), is.finite(value))

scale_settings <- tibble(
  variable = c("pVL", "DNA", "CD4", "CD8"),
  value_scale_used = c(
    "raw copies/mL -> log10",
    "already log10 copies per 10^6 CD4 cells",
    "raw cells/mL -> log10",
    "raw cells/mL -> log10"
  ),
  lod_scale_used = c(
    "raw copies/mL -> log10",
    "no LOD used",
    "no LOD used",
    "no LOD used"
  )
)

print(scale_settings)

write.csv(
  scale_settings,
  file.path(out_dir, "scale_settings_used.csv"),
  row.names = FALSE
)

scale_check <- data_fit %>%
  group_by(variable) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    n = n(),
    n_censored = sum(is_censored, na.rm = TRUE),
    .groups = "drop"
  )

print(scale_check)

write.csv(
  scale_check,
  file.path(out_dir, "data_scale_check.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Residual scaling
# ============================================================

if (USE_VARIABLE_SCALING) {
  residual_scale_table <- data_fit %>%
    group_by(variable) %>%
    summarise(
      residual_scale = sd(value[!is_censored], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      residual_scale = ifelse(
        !is.finite(residual_scale) | residual_scale <= 0,
        1,
        residual_scale
      )
    )
} else {
  residual_scale_table <- tibble(
    variable = FIT_VARIABLES,
    residual_scale = 1
  )
}

print(residual_scale_table)

write.csv(
  residual_scale_table,
  file.path(out_dir, "residual_scale_table.csv"),
  row.names = FALSE
)

# ============================================================
# 9. ODE model
# ============================================================

ati_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    T <- max(T, EPS)
    L <- max(L, EPS)
    I <- max(I, EPS)
    V <- max(V, EPS)
    C <- max(C, EPS)
    
    K_C <- max(K_C, EPS)
    K_kill <- max(K_kill, EPS)
    xi_C <- max(xi_C, EPS)
    
    C_scaled <- C / CD8_scale
    
    CD8_activation_fraction <- I / (K_C + I)
    CD8_killing_saturation <- C_scaled / (K_kill + C_scaled)
    CD8_killing_rate <- k_C_G * CD8_killing_saturation
    
    # xi_C retained:
    # This is the non-cytolytic CD8 suppression on viral production.
    p_eff <- p_G / (1 + xi_C * C_scaled)
    
    dT <- d_T * T0 -
      d_T * T -
      beta * V * T
    
    dL <- f_L * beta * V * T -
      d_L * L -
      a_L * L
    
    dI <- (1 - f_L) * beta * V * T +
      a_L * L -
      delta0 * I -
      CD8_killing_rate * I
    
    dV <- p_eff * I -
      c * V
    
    dC <- d_C * (C0 - C) +
      rho_C * CD8_activation_fraction * C
    
    list(c(dT, dL, dI, dV, dC))
  })
}

# ============================================================
# 10. Simulation and prediction functions
# ============================================================

simulate_one <- function(id_row, pars, times) {
  
  T0 <- max(as.numeric(id_row$T0), EPS)
  L0 <- max(as.numeric(id_row$L0), EPS)
  I0_init <- max(as.numeric(id_row$I0), EPS)
  V0 <- max(as.numeric(id_row$V0), EPS)
  C0 <- max(as.numeric(id_row$C0), EPS)
  
  p_G <- as.numeric(pars["p"])
  k_C_G <- as.numeric(pars["k_C"])
  
  I0 <- I0_init
  
  if (RECOMPUTE_I0_FROM_QSS) {
    I0 <- max(as.numeric(pars["c"]) * V0 / max(p_G, EPS), EPS)
  }
  
  state <- c(
    T = T0,
    L = L0,
    I = I0,
    V = V0,
    C = C0
  )
  
  times <- sort(unique(c(0, times)))
  
  out <- tryCatch({
    ode(
      y = state,
      times = times,
      func = ati_model,
      parms = c(
        pars,
        T0 = T0,
        C0 = C0,
        p_G = p_G,
        k_C_G = k_C_G
      ),
      method = "lsoda"
    )
  }, error = function(e) {
    message("ODE error for id ", id_row$id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) {
    return(NULL)
  }
  
  sim <- as.data.frame(out)
  
  sim %>%
    mutate(
      id = id_row$id,
      group = id_row$group,
      G_early = ifelse(is_early_group(id_row$group), 1, 0),
      T0 = T0,
      L0 = L0,
      I0 = I0,
      I0_init_table = I0_init,
      V0 = V0,
      C0 = C0,
      C0_source = id_row$C0_source,
      a_L = as.numeric(pars["a_L"]),
      p_G = p_G,
      k_C_G = k_C_G,
      xi_C = as.numeric(pars["xi_C"]),
      rho_C = as.numeric(pars["rho_C"]),
      delta0 = as.numeric(pars["delta0"]),
      virion_clearance_c = as.numeric(pars["c"]),
      d_C_value = as.numeric(pars["d_C"]),
      K_C_value = as.numeric(pars["K_C"]),
      K_kill_value = as.numeric(pars["K_kill"]),
      C_scaled = C / as.numeric(pars["CD8_scale"]),
      CD8_activation_fraction = I / (as.numeric(pars["K_C"]) + I),
      CD8_killing_saturation = C_scaled / (as.numeric(pars["K_kill"]) + C_scaled),
      CD8_killing_rate = k_C_G * CD8_killing_saturation,
      p_eff = p_G / (1 + as.numeric(pars["xi_C"]) * C_scaled),
      p_eff_over_p_G = p_eff / p_G,
      CD4_total = T + L + I,
      CD8_total = C,
      pred_pVL_log10 = safe_log10(V),
      pred_DNA_log10 = safe_log10(1e6 * (L + I) / pmax(CD4_total, EPS)),
      pred_CD4_log10 = safe_log10(CD4_total),
      pred_CD8_log10 = safe_log10(CD8_total)
    )
}

simulate_all <- function(pars) {
  
  sim_list <- vector("list", nrow(init_fit))
  
  for (i in seq_len(nrow(init_fit))) {
    
    id_i <- init_fit$id[i]
    
    times_i <- data_fit %>%
      filter(id == id_i, !is.na(time)) %>%
      pull(time) %>%
      unique()
    
    times_i <- sort(unique(c(0, times_i)))
    
    sim_i <- simulate_one(
      id_row = init_fit[i, ],
      pars = pars,
      times = times_i
    )
    
    if (is.null(sim_i)) {
      return(NULL)
    }
    
    sim_list[[i]] <- sim_i
  }
  
  bind_rows(sim_list)
}

make_predictions <- function(sim) {
  
  pred_cols <- c(
    pred_pVL_log10 = "pVL",
    pred_DNA_log10 = "DNA",
    pred_CD4_log10 = "CD4",
    pred_CD8_log10 = "CD8"
  )
  
  selected_pred_cols <- names(pred_cols)[pred_cols %in% FIT_VARIABLES]
  
  sim %>%
    select(id, group, time, all_of(selected_pred_cols)) %>%
    pivot_longer(
      cols = all_of(selected_pred_cols),
      names_to = "variable",
      values_to = "pred"
    ) %>%
    mutate(
      variable = recode(
        variable,
        pred_pVL_log10 = "pVL",
        pred_DNA_log10 = "DNA",
        pred_CD4_log10 = "CD4",
        pred_CD8_log10 = "CD8"
      )
    )
}

# ============================================================
# 11. Fitting functions
# ============================================================

default_start_free_pars <- c(
  a_L = 1e-3,
  p = 6e3,
  k_C = 0.01,
  xi_C = 0.5,
  rho_C = 0.05
)

make_start_vector <- function(start_free_pars = default_start_free_pars) {
  c(
    log_a_L = log(unname(start_free_pars["a_L"])),
    log_p = log(unname(start_free_pars["p"])),
    log_k_C = log(unname(start_free_pars["k_C"])),
    log_xi_C = log(unname(start_free_pars["xi_C"])),
    log_rho_C = log(unname(start_free_pars["rho_C"]))
  )
}

make_lower_vector <- function() {
  c(
    log_a_L = log(1e-8),
    log_p = log(1),
    log_k_C = log(1e-6),
    log_xi_C = log(1e-5),
    log_rho_C = log(1e-6)
  )
}

make_upper_vector <- function() {
  c(
    log_a_L = log(1),
    log_p = log(1e8),
    log_k_C = log(5),
    log_xi_C = log(100),
    log_rho_C = log(1)
  )
}

unpack_pars <- function(x, fixed_pars) {
  c(
    fixed_pars,
    a_L = exp(unname(x["log_a_L"])),
    p = exp(unname(x["log_p"])),
    k_C = exp(unname(x["log_k_C"])),
    xi_C = exp(unname(x["log_xi_C"])),
    rho_C = exp(unname(x["log_rho_C"]))
  )
}

make_residual_table <- function(x, fixed_pars) {
  
  pars <- unpack_pars(x, fixed_pars)
  
  sim <- simulate_all(pars = pars)
  
  if (is.null(sim)) {
    return(NULL)
  }
  
  pred <- make_predictions(sim)
  
  joined <- data_fit %>%
    left_join(pred, by = c("id", "group", "time", "variable")) %>%
    left_join(residual_scale_table, by = "variable")
  
  if (
    any(is.na(joined$pred)) ||
    any(!is.finite(joined$pred)) ||
    any(is.na(joined$residual_scale)) ||
    any(!is.finite(joined$residual_scale))
  ) {
    return(NULL)
  }
  
  joined %>%
    mutate(
      raw_resid = case_when(
        is_censored & pred <= lod ~ 0,
        is_censored & pred > lod ~ pred - lod,
        TRUE ~ pred - value
      ),
      resid = raw_resid / residual_scale
    )
}

make_residuals <- function(x, fixed_pars) {
  
  joined <- make_residual_table(x = x, fixed_pars = fixed_pars)
  
  if (is.null(joined)) {
    return(rep(1e6, nrow(data_fit)))
  }
  
  resid <- joined$resid
  
  if (any(!is.finite(resid))) {
    return(rep(1e6, nrow(data_fit)))
  }
  
  resid
}

make_random_start <- function() {
  c(
    a_L = exp(runif(1, log(1e-5), log(1e-1))),
    p = exp(runif(1, log(1e2), log(1e5))),
    k_C = exp(runif(1, log(1e-4), log(1))),
    xi_C = exp(runif(1, log(1e-3), log(10))),
    rho_C = exp(runif(1, log(1e-4), log(0.5)))
  )
}

make_random_start_set <- function(n_start = 5, seed = 123) {
  
  set.seed(seed)
  
  start_set <- list()
  start_set[[length(start_set) + 1]] <- default_start_free_pars
  
  for (i in seq_len(n_start)) {
    start_set[[length(start_set) + 1]] <- make_random_start()
  }
  
  start_set
}

fit_one_parameter_setting <- function(parameter_row, start_free_pars, start_id = NA_integer_) {
  
  fixed_pars <- row_to_fixed_pars(parameter_row)
  
  x0 <- make_start_vector(start_free_pars = start_free_pars)
  
  fit <- tryCatch({
    nls.lm(
      par = x0,
      lower = make_lower_vector(),
      upper = make_upper_vector(),
      fn = make_residuals,
      fixed_pars = fixed_pars,
      control = nls.lm.control(
        maxiter = 300,
        ftol = 1e-10,
        ptol = 1e-10,
        gtol = 1e-10
      )
    )
  }, error = function(e) {
    message("Fitting failed for setting ", parameter_row$setting_id, ": ", e$message)
    return(NULL)
  })
  
  failed_summary <- function(convergence_info = NA_integer_) {
    tibble(
      setting_id = parameter_row$setting_id,
      description = parameter_row$description,
      n = nrow(data_fit),
      q = length(x0),
      RSS = NA_real_,
      raw_RSS = NA_real_,
      convergence_info = convergence_info,
      start_id = start_id
    )
  }
  
  failed_params <- function() {
    tibble(
      setting_id = parameter_row$setting_id,
      description = parameter_row$description,
      a_L = NA_real_,
      p = NA_real_,
      k_C = NA_real_,
      xi_C = NA_real_,
      rho_C = NA_real_,
      start_id = start_id
    )
  }
  
  if (is.null(fit)) {
    return(list(fit = NULL, summary = failed_summary(), params = failed_params()))
  }
  
  x_hat <- fit$par
  pars_hat <- unpack_pars(x_hat, fixed_pars)
  
  residual_table_hat <- make_residual_table(x = x_hat, fixed_pars = fixed_pars)
  
  if (is.null(residual_table_hat)) {
    return(list(
      fit = fit,
      summary = failed_summary(convergence_info = fit$info),
      params = failed_params()
    ))
  }
  
  resid_hat <- residual_table_hat$resid
  raw_resid_hat <- residual_table_hat$raw_resid
  
  RSS <- sum(resid_hat^2)
  raw_RSS <- sum(raw_resid_hat^2)
  
  summary_row <- tibble(
    setting_id = parameter_row$setting_id,
    description = parameter_row$description,
    n = length(resid_hat),
    q = length(x_hat),
    RSS = RSS,
    raw_RSS = raw_RSS,
    convergence_info = fit$info,
    start_id = start_id
  )
  
  param_row <- tibble(
    setting_id = parameter_row$setting_id,
    description = parameter_row$description,
    a_L = unname(pars_hat["a_L"]),
    p = unname(pars_hat["p"]),
    k_C = unname(pars_hat["k_C"]),
    xi_C = unname(pars_hat["xi_C"]),
    rho_C = unname(pars_hat["rho_C"]),
    start_id = start_id
  )
  
  list(
    fit = fit,
    summary = summary_row,
    params = param_row
  )
}

fit_one_parameter_setting_multistart <- function(parameter_row, n_start = 5, seed = 123) {
  
  start_set <- make_random_start_set(n_start = n_start, seed = seed)
  
  fit_candidates <- purrr::map(
    seq_along(start_set),
    function(i) {
      cat(
        "\nFixed-parameter setting:", parameter_row$setting_id,
        "|", parameter_row$description,
        "| start:", i,
        "\n"
      )
      
      fit_one_parameter_setting(
        parameter_row = parameter_row,
        start_free_pars = start_set[[i]],
        start_id = i
      )
    }
  )
  
  summary_candidates <- purrr::map_dfr(fit_candidates, "summary") %>%
    mutate(n_random_start = n_start)
  
  valid_candidates <- summary_candidates %>%
    filter(!is.na(RSS), is.finite(RSS))
  
  if (nrow(valid_candidates) == 0) {
    failed_summary <- summary_candidates %>%
      slice(1) %>%
      mutate(best_start_id = NA_integer_)
    
    failed_params <- purrr::map_dfr(fit_candidates, "params") %>%
      slice(1) %>%
      mutate(best_start_id = NA_integer_, n_random_start = n_start)
    
    return(list(
      best_fit = list(fit = NULL, summary = failed_summary, params = failed_params),
      all_summaries = summary_candidates
    ))
  }
  
  best_start_id <- valid_candidates %>%
    arrange(RSS) %>%
    slice(1) %>%
    pull(start_id)
  
  best_fit <- fit_candidates[[best_start_id]]
  
  best_fit$summary <- best_fit$summary %>%
    mutate(best_start_id = best_start_id, n_random_start = n_start)
  
  best_fit$params <- best_fit$params %>%
    mutate(best_start_id = best_start_id, n_random_start = n_start)
  
  list(
    best_fit = best_fit,
    all_summaries = summary_candidates
  )
}

# ============================================================
# 12. Run fixed-parameter calibration
# ============================================================

calibration_fit_list <- vector("list", nrow(fixed_parameter_grid))

for (i in seq_len(nrow(fixed_parameter_grid))) {
  calibration_fit_list[[i]] <- fit_one_parameter_setting_multistart(
    parameter_row = fixed_parameter_grid[i, ],
    n_start = N_RANDOM_START,
    seed = RANDOM_SEED_BASE + i
  )
}

fit_list_best <- purrr::map(calibration_fit_list, "best_fit")

parameter_fit_table <- purrr::map_dfr(fit_list_best, "summary")

best_fit_params <- purrr::map_dfr(fit_list_best, "params")

all_random_fit_summaries <- purrr::map_dfr(
  seq_along(calibration_fit_list),
  function(i) {
    calibration_fit_list[[i]]$all_summaries
  }
)

write.csv(
  all_random_fit_summaries,
  file.path(out_dir, "all_random_fit_summaries.csv"),
  row.names = FALSE
)

# ============================================================
# 13. Diagnostics for each fixed-parameter setting
# ============================================================

setting_diagnostics_list <- vector("list", nrow(fixed_parameter_grid))
v_production_summary_list <- vector("list", nrow(fixed_parameter_grid))

for (i in seq_len(nrow(fixed_parameter_grid))) {
  
  best_fit_i <- fit_list_best[[i]]
  parameter_row_i <- fixed_parameter_grid[i, ]
  fixed_pars_i <- row_to_fixed_pars(parameter_row_i)
  
  if (is.null(best_fit_i$fit)) {
    setting_diagnostics_list[[i]] <- tibble(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description,
      variable = FIT_VARIABLES,
      n = NA_integer_,
      n_censored = NA_integer_,
      mean_raw_resid = NA_real_,
      median_raw_resid = NA_real_,
      RMSE_raw = NA_real_,
      weighted_RSS = NA_real_
    )
    
    v_production_summary_list[[i]] <- tibble(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description,
      group = NA_character_,
      p = NA_real_,
      xi_C = NA_real_,
      median_C_scaled = NA_real_,
      median_p_eff = NA_real_,
      median_p_eff_over_p_G = NA_real_,
      median_pVL_log10 = NA_real_
    )
    next
  }
  
  residual_table_i <- make_residual_table(
    x = best_fit_i$fit$par,
    fixed_pars = fixed_pars_i
  )
  
  pars_i <- unpack_pars(best_fit_i$fit$par, fixed_pars_i)
  sim_i <- simulate_all(pars = pars_i)
  
  setting_diagnostics_list[[i]] <- residual_table_i %>%
    group_by(variable) %>%
    summarise(
      n = n(),
      n_censored = sum(is_censored, na.rm = TRUE),
      mean_raw_resid = mean(raw_resid, na.rm = TRUE),
      median_raw_resid = median(raw_resid, na.rm = TRUE),
      RMSE_raw = sqrt(mean(raw_resid^2, na.rm = TRUE)),
      weighted_RSS = sum(resid^2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description
    ) %>%
    select(setting_id, description, variable, everything())
  
  v_production_summary_list[[i]] <- sim_i %>%
    group_by(group) %>%
    summarise(
      p = median(p_G, na.rm = TRUE),
      xi_C = median(xi_C, na.rm = TRUE),
      median_C_scaled = median(C_scaled, na.rm = TRUE),
      median_p_eff = median(p_eff, na.rm = TRUE),
      median_p_eff_over_p_G = median(p_eff_over_p_G, na.rm = TRUE),
      median_pVL_log10 = median(pred_pVL_log10, na.rm = TRUE),
      median_CD8_killing_rate = median(CD8_killing_rate, na.rm = TRUE),
      median_CD8_activation_fraction = median(CD8_activation_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description
    ) %>%
    select(setting_id, description, group, everything())
}

setting_diagnostics_by_variable <- bind_rows(setting_diagnostics_list)
v_production_summary_by_setting <- bind_rows(v_production_summary_list)

write.csv(
  setting_diagnostics_by_variable,
  file.path(out_dir, "fixed_parameter_diagnostics_by_variable.csv"),
  row.names = FALSE
)

write.csv(
  v_production_summary_by_setting,
  file.path(out_dir, "v_production_effect_summary_by_setting.csv"),
  row.names = FALSE
)

setting_bias_wide <- setting_diagnostics_by_variable %>%
  select(setting_id, variable, mean_raw_resid, RMSE_raw, weighted_RSS) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(mean_raw_resid, RMSE_raw, weighted_RSS),
    names_sep = "_"
  )

v_production_wide <- v_production_summary_by_setting %>%
  group_by(setting_id) %>%
  summarise(
    median_p = median(p, na.rm = TRUE),
    median_xi_C = median(xi_C, na.rm = TRUE),
    median_p_eff = median(median_p_eff, na.rm = TRUE),
    median_p_eff_over_p_G = median(median_p_eff_over_p_G, na.rm = TRUE),
    median_CD8_killing_rate = median(median_CD8_killing_rate, na.rm = TRUE),
    median_CD8_activation_fraction = median(median_CD8_activation_fraction, na.rm = TRUE),
    .groups = "drop"
  )

fixed_parameter_ranking <- parameter_fit_table %>%
  left_join(fixed_parameter_grid, by = c("setting_id", "description")) %>%
  left_join(best_fit_params, by = c("setting_id", "description", "start_id", "best_start_id", "n_random_start")) %>%
  left_join(setting_bias_wide, by = "setting_id") %>%
  left_join(v_production_wide, by = "setting_id") %>%
  mutate(
    mean_raw_resid_pVL_for_penalty = replace_na(mean_raw_resid_pVL, 0),
    mean_raw_resid_CD8_for_penalty = replace_na(mean_raw_resid_CD8, 0),
    pVL_underprediction_penalty = pmax(0, -mean_raw_resid_pVL_for_penalty)^2,
    CD8_overprediction_penalty = pmax(0, mean_raw_resid_CD8_for_penalty)^2,
    known_problem_penalty = 20 * pVL_underprediction_penalty + 20 * CD8_overprediction_penalty,
    problem_score = RSS + known_problem_penalty,
    rank_by_RSS = rank(RSS, ties.method = "first", na.last = "keep"),
    rank_by_problem_score = rank(problem_score, ties.method = "first", na.last = "keep")
  ) %>%
  arrange(problem_score, RSS)

print(fixed_parameter_ranking)

write.csv(
  fixed_parameter_ranking,
  file.path(out_dir, "fixed_parameter_ranking.csv"),
  row.names = FALSE
)

write.csv(
  best_fit_params,
  file.path(out_dir, "fixed_parameter_best_fit_free_params.csv"),
  row.names = FALSE
)

write.csv(
  parameter_fit_table,
  file.path(out_dir, "fixed_parameter_fit_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 14. Extract best setting and save detailed outputs
# ============================================================

if (all(is.na(fixed_parameter_ranking$problem_score))) {
  stop("All parameter-setting fits failed. Check all_random_fit_summaries.csv and input scale checks.")
}

best_setting_id <- fixed_parameter_ranking %>%
  arrange(problem_score, RSS) %>%
  slice(1) %>%
  pull(setting_id)

best_setting_index <- match(best_setting_id, fixed_parameter_grid$setting_id)
best_parameter_row <- fixed_parameter_grid[best_setting_index, ]
best_fixed_pars <- row_to_fixed_pars(best_parameter_row)
best_fit_overall <- fit_list_best[[best_setting_index]]

best_residual_table <- make_residual_table(
  x = best_fit_overall$fit$par,
  fixed_pars = best_fixed_pars
)

best_pars <- unpack_pars(
  x = best_fit_overall$fit$par,
  fixed_pars = best_fixed_pars
)

best_sim <- simulate_all(pars = best_pars)
best_pred <- make_predictions(best_sim)

write.csv(
  best_residual_table,
  file.path(out_dir, "best_setting_residual_table.csv"),
  row.names = FALSE
)

write.csv(
  best_sim,
  file.path(out_dir, "best_setting_simulation_full.csv"),
  row.names = FALSE
)

write.csv(
  best_pred,
  file.path(out_dir, "best_setting_predictions_long.csv"),
  row.names = FALSE
)

best_setting_params <- tibble(
  setting_id = best_parameter_row$setting_id,
  description = best_parameter_row$description,
  d_T = best_pars["d_T"],
  beta = best_pars["beta"],
  f_L = best_pars["f_L"],
  d_L = best_pars["d_L"],
  delta0 = best_pars["delta0"],
  c = best_pars["c"],
  d_C = best_pars["d_C"],
  K_C = best_pars["K_C"],
  K_kill = best_pars["K_kill"],
  CD8_scale = best_pars["CD8_scale"],
  a_L = best_pars["a_L"],
  p = best_pars["p"],
  k_C = best_pars["k_C"],
  xi_C = best_pars["xi_C"],
  rho_C = best_pars["rho_C"],
  recompute_I0_from_qss = RECOMPUTE_I0_FROM_QSS
)

print(best_setting_params)

write.csv(
  best_setting_params,
  file.path(out_dir, "best_setting_parameters_to_use_before_model_selection.csv"),
  row.names = FALSE
)

best_rss_by_variable <- best_residual_table %>%
  group_by(variable) %>%
  summarise(
    n = n(),
    n_censored = sum(is_censored, na.rm = TRUE),
    mean_raw_resid = mean(raw_resid, na.rm = TRUE),
    median_raw_resid = median(raw_resid, na.rm = TRUE),
    RMSE_raw = sqrt(mean(raw_resid^2, na.rm = TRUE)),
    weighted_RSS = sum(resid^2, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  best_rss_by_variable,
  file.path(out_dir, "best_setting_RSS_by_variable.csv"),
  row.names = FALSE
)

best_group_bias_summary <- best_residual_table %>%
  group_by(variable, group) %>%
  summarise(
    n = n(),
    n_censored = sum(is_censored, na.rm = TRUE),
    mean_raw_resid = mean(raw_resid, na.rm = TRUE),
    median_raw_resid = median(raw_resid, na.rm = TRUE),
    RMSE_raw = sqrt(mean(raw_resid^2, na.rm = TRUE)),
    mean_weighted_resid = mean(resid, na.rm = TRUE),
    median_weighted_resid = median(resid, na.rm = TRUE),
    RMSE_weighted = sqrt(mean(resid^2, na.rm = TRUE)),
    weighted_RSS = sum(resid^2, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  best_group_bias_summary,
  file.path(out_dir, "best_setting_group_bias_summary.csv"),
  row.names = FALSE
)

best_cd8_effect_summary <- best_sim %>%
  group_by(id, group) %>%
  summarise(
    C0 = first(C0),
    C0_source = first(C0_source),
    median_C = median(C, na.rm = TRUE),
    median_C_scaled = median(C_scaled, na.rm = TRUE),
    median_CD8_activation_fraction = median(CD8_activation_fraction, na.rm = TRUE),
    median_CD8_killing_saturation = median(CD8_killing_saturation, na.rm = TRUE),
    median_CD8_killing_rate = median(CD8_killing_rate, na.rm = TRUE),
    median_p_eff = median(p_eff, na.rm = TRUE),
    median_p_eff_over_p_G = median(p_eff_over_p_G, na.rm = TRUE),
    a_L = median(a_L, na.rm = TRUE),
    p_G = median(p_G, na.rm = TRUE),
    k_C_G = median(k_C_G, na.rm = TRUE),
    xi_C = median(xi_C, na.rm = TRUE),
    rho_C = median(rho_C, na.rm = TRUE),
    delta0 = median(delta0, na.rm = TRUE),
    c = median(virion_clearance_c, na.rm = TRUE),
    d_C = median(d_C_value, na.rm = TRUE),
    K_C = median(K_C_value, na.rm = TRUE),
    K_kill = median(K_kill_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(group, id)

write.csv(
  best_cd8_effect_summary,
  file.path(out_dir, "best_setting_CD8_effect_summary_by_id.csv"),
  row.names = FALSE
)

# ============================================================
# 15. Figures
# ============================================================

p_parameter_ranking <- fixed_parameter_ranking %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  mutate(setting_label = factor(setting_label, levels = rev(setting_label))) %>%
  ggplot(aes(x = setting_label, y = problem_score)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Fixed-parameter calibration ranking with xi_C",
    x = "Fixed-parameter setting",
    y = "Problem score = weighted RSS + pVL/CD8 bias penalty"
  )

print(p_parameter_ranking)
save_plot_png_pdf(p_parameter_ranking, "Fig1_fixed_parameter_ranking_with_xiC", width = 10, height = 7)

p_obs_pred <- best_residual_table %>%
  ggplot(aes(x = value, y = pred)) +
  geom_point(aes(shape = is_censored), alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  labs(
    title = paste0("Observed vs predicted: best fixed setting ", best_setting_id),
    x = "Observed value used for fitting",
    y = "Predicted value",
    shape = "Censored"
  )

print(p_obs_pred)
save_plot_png_pdf(p_obs_pred, "Fig2_observed_vs_predicted_best_setting", width = 9, height = 5)

p_resid <- best_residual_table %>%
  ggplot(aes(x = variable, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +
  theme_bw() +
  labs(
    title = paste0("Weighted residuals by variable: best fixed setting ", best_setting_id),
    x = "Variable",
    y = "Weighted residual"
  )

print(p_resid)
save_plot_png_pdf(p_resid, "Fig3_residual_by_variable_best_setting", width = 7, height = 5)

p_timecourse <- best_residual_table %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = pred, group = id), alpha = 0.35, linewidth = 0.5) +
  geom_point(aes(y = value, shape = is_censored), alpha = 0.7, size = 1.8) +
  facet_grid(variable ~ group, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste0("Observed and predicted time courses: best fixed setting ", best_setting_id),
    x = "Time",
    y = "Value used for fitting",
    shape = "Censored"
  )

print(p_timecourse)
save_plot_png_pdf(p_timecourse, "Fig4_timecourse_observed_predicted_best_setting", width = 10, height = 8)

p_cd8_time <- best_sim %>%
  select(
    id,
    group,
    time,
    C,
    C_scaled,
    CD8_activation_fraction,
    CD8_killing_saturation,
    CD8_killing_rate,
    p_eff_over_p_G
  ) %>%
  pivot_longer(
    cols = c(C, C_scaled, CD8_activation_fraction, CD8_killing_saturation, CD8_killing_rate, p_eff_over_p_G),
    names_to = "quantity",
    values_to = "value"
  ) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  facet_grid(quantity ~ group, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste0("CD8 state, killing, and p suppression: best fixed setting ", best_setting_id),
    x = "Time",
    y = "Value"
  )

print(p_cd8_time)
save_plot_png_pdf(p_cd8_time, "Fig5_CD8_states_killing_and_xiC_suppression_best_setting", width = 11, height = 9)

p_bias_compare <- setting_diagnostics_by_variable %>%
  filter(variable %in% c("pVL", "CD8")) %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  ggplot(aes(x = setting_label, y = mean_raw_resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Mean raw residual bias for pVL and CD8 across fixed-parameter settings",
    x = "Fixed-parameter setting",
    y = "Mean raw residual: predicted - observed"
  )

print(p_bias_compare)
save_plot_png_pdf(p_bias_compare, "Fig6_pVL_CD8_bias_across_fixed_settings", width = 11, height = 7)

p_xiC_compare <- v_production_summary_by_setting %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  ggplot(aes(x = setting_label, y = median_p_eff_over_p_G)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ group) +
  theme_bw() +
  labs(
    title = "Effective viral production fraction across fixed-parameter settings",
    x = "Fixed-parameter setting",
    y = "Median p_eff / p"
  )

print(p_xiC_compare)
save_plot_png_pdf(p_xiC_compare, "Fig7_p_eff_over_p_across_fixed_settings", width = 11, height = 7)

# ============================================================
# 16. Save run settings and final message
# ============================================================

run_settings <- tibble(
  item = c(
    "N_RANDOM_START",
    "RANDOM_SEED_BASE",
    "USE_VARIABLE_SCALING",
    "FIT_VARIABLES",
    "RECOMPUTE_I0_FROM_QSS",
    "data_file",
    "init_file",
    "CD8_scale_reference",
    "fitted_common_parameters",
    "tested_fixed_parameters",
    "xi_C_structure",
    "ranking_rule",
    "out_dir"
  ),
  value = c(
    as.character(N_RANDOM_START),
    as.character(RANDOM_SEED_BASE),
    as.character(USE_VARIABLE_SCALING),
    paste(FIT_VARIABLES, collapse = ","),
    as.character(RECOMPUTE_I0_FROM_QSS),
    DATA_FILE,
    INIT_FILE,
    as.character(CD8_SCALE),
    "a_L, p, k_C, xi_C, rho_C",
    "delta0, c, d_C, K_C, K_kill",
    "p_eff = p / (1 + xi_C * C_scaled), C_scaled = C / 1e6",
    "arrange by problem_score first, then RSS; problem_score penalizes pVL underprediction and CD8 overprediction",
    out_dir
  )
)

write.csv(
  run_settings,
  file.path(out_dir, "run_settings.csv"),
  row.names = FALSE
)

cat("\nFixed-parameter calibration with xi_C finished.\n")
cat("\nBest fixed-parameter setting by problem_score:\n")
cat(best_setting_id, " - ", best_parameter_row$description, "\n", sep = "")

cat("\nBest parameters to consider before model selection:\n")
print(best_setting_params)

cat("\nMain files to check:\n")
cat("1. fixed_parameter_test_grid.csv\n")
cat("2. fixed_parameter_ranking.csv\n")
cat("3. fixed_parameter_diagnostics_by_variable.csv\n")
cat("4. v_production_effect_summary_by_setting.csv\n")
cat("5. fixed_parameter_best_fit_free_params.csv\n")
cat("6. best_setting_parameters_to_use_before_model_selection.csv\n")
cat("7. best_setting_RSS_by_variable.csv\n")
cat("8. best_setting_group_bias_summary.csv\n")
cat("9. best_setting_CD8_effect_summary_by_id.csv\n")
cat("10. best_setting_residual_table.csv\n")
cat("11. best_setting_simulation_full.csv\n")
cat("12. best_setting_predictions_long.csv\n")

cat("\nMain figures to check:\n")
cat("1. Fig1_fixed_parameter_ranking_with_xiC.png\n")
cat("2. Fig2_observed_vs_predicted_best_setting.png\n")
cat("3. Fig3_residual_by_variable_best_setting.png\n")
cat("4. Fig4_timecourse_observed_predicted_best_setting.png\n")
cat("5. Fig5_CD8_states_killing_and_xiC_suppression_best_setting.png\n")
cat("6. Fig6_pVL_CD8_bias_across_fixed_settings.png\n")
cat("7. Fig7_p_eff_over_p_across_fixed_settings.png\n")

cat("\nAll results saved to:\n")
cat(out_dir, "\n")

cat("\nFigures saved to:\n")
cat(fig_dir, "\n")
