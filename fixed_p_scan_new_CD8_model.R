library(deSolve)
library(tidyverse)
library(minpack.lm)

# ============================================================
# Fixed-p scan for the new CD8 model
# ============================================================
# Purpose:
#   Fix viral production rate p at several candidate values.
#   For each fixed p value, fit a_L, k_C, xi_C, and rho_C.
#   This is NOT model selection.
#
# New CD8 equation:
#   dC/dt = d_C * (C_base_group - C) + rho_C * I/(K_C + I) * C
#
# Viral production with non-cytolytic CD8 suppression:
#   p_eff = p / (1 + xi_C * C_scaled)
#   C_scaled = C / 1e6
#
# For each fixed p value, this script fits:
#   a_L, k_C, xi_C, rho_C
#
# Fixed biological background used here:
#   S03-like setting: delta0 = 0.45, c = 23, d_C = 0.75,
#   K_C = 1e4, K_kill = 5.
# ============================================================

# ============================================================
# 0. Global settings
# ============================================================

EPS <- 1e-12

N_RANDOM_START <- 10
RANDOM_SEED_BASE <- 1000

USE_VARIABLE_SCALING <- TRUE
FIT_VARIABLES <- c("pVL", "DNA", "CD4", "CD8")

# Keep FALSE for the first comparison.
# If TRUE, I0 is recalculated as c*V0/p for each fitted p and fixed c.
RECOMPUTE_I0_FROM_QSS <- FALSE

CD8_SCALE <- 1e6

# Use your usual project data folder.
# If your files are in the same folder as the script, change this to getwd().
DATA_DIR <- "D:/R-4.6.0/R-ATI-project/Rdata"

DATA_FILE <- file.path(DATA_DIR, "raw_data_with_cd8_counts.csv")
INIT_FILE <- file.path(DATA_DIR, "init_table_total_CD8_no_exhaustion.csv")

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  DATA_DIR,
  "fixed_p_scan_new_CD8_model",
  run_tag
)

fig_dir <- file.path(out_dir, "figures")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nDATA_FILE:\n", DATA_FILE, "\n", sep = "")
cat("DATA_FILE exists: ", file.exists(DATA_FILE), "\n", sep = "")
cat("\nINIT_FILE:\n", INIT_FILE, "\n", sep = "")
cat("INIT_FILE exists: ", file.exists(INIT_FILE), "\n", sep = "")

if (!file.exists(DATA_FILE)) {
  stop("找不到 DATA_FILE: ", DATA_FILE, "\n当前工作目录是: ", getwd())
}

if (!file.exists(INIT_FILE)) {
  stop("找不到 INIT_FILE: ", INIT_FILE, "\n当前工作目录是: ", getwd())
}

# ============================================================
# 1. Selected fixed-parameter sets
# ============================================================
# S03: biologically safer mainline candidate, keeps c = 23 and uses moderate K values.
# S05: best previous numerical fit, but c = 3 and fitted p may become very low.
# S10: second previous numerical candidate, keeps c = 23 and uses slower CD8 homeostatic return.

fixed_parameter_grid <- tribble(
  ~setting_id, ~description,         ~p_fixed, ~d_T,  ~beta,     ~f_L,  ~d_L,  ~delta0, ~c,  ~d_C, ~K_C,  ~K_kill, ~CD8_scale,
  "P4000",    "fixed_p_4000",       4000,     0.01, 1.58e-8,   0.005, 5e-4,  0.45,    23,  0.75, 1e4,   5.0,     CD8_SCALE,
  "P6000",    "fixed_p_6000_main",  6000,     0.01, 1.58e-8,   0.005, 5e-4,  0.45,    23,  0.75, 1e4,   5.0,     CD8_SCALE,
  "P8000",    "fixed_p_8000",       8000,     0.01, 1.58e-8,   0.005, 5e-4,  0.45,    23,  0.75, 1e4,   5.0,     CD8_SCALE
)

write.csv(
  fixed_parameter_grid,
  file.path(out_dir, "fixed_p_test_grid.csv"),
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
    p = as.numeric(parameter_row$p_fixed),
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
# 6. Add C0 and group-level CD8 baseline
# ============================================================
# C0 is the individual initial CD8 value for the simulation state.
# C_base_group is the group-level homeostatic baseline used in the CD8 equation.

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
    ),
    C_base_group = case_when(
      is.finite(CD8_count_group) ~ CD8_count_group,
      TRUE ~ overall_cd8_median
    ),
    C_base_source = case_when(
      is.finite(CD8_count_group) ~ "group_median_CD8_all_times",
      TRUE ~ "overall_median_CD8_all_times"
    )
  ) %>%
  filter(
    is.finite(C0),
    C0 > 0,
    is.finite(C_base_group),
    C_base_group > 0
  )

write.csv(
  init_fit,
  file.path(out_dir, "init_table_with_C0_and_Cbase_used.csv"),
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

    # New CD8 equation:
    # C0 is only the initial value.
    # C_base_group is the homeostatic set point.
    dC <- d_C * (C_base_group - C) +
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
  C_base_group <- max(as.numeric(id_row$C_base_group), EPS)

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
        C_base_group = C_base_group,
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
      C_base_group = C_base_group,
      C_base_source = id_row$C_base_source,
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
  k_C = 0.01,
  xi_C = 0.5,
  rho_C = 0.05
)

make_start_vector <- function(start_free_pars = default_start_free_pars) {
  c(
    log_a_L = log(unname(start_free_pars["a_L"])),
    log_k_C = log(unname(start_free_pars["k_C"])),
    log_xi_C = log(unname(start_free_pars["xi_C"])),
    log_rho_C = log(unname(start_free_pars["rho_C"]))
  )
}

make_lower_vector <- function() {
  c(
    log_a_L = log(1e-8),
    log_k_C = log(1e-6),
    log_xi_C = log(1e-5),
    log_rho_C = log(1e-6)
  )
}

make_upper_vector <- function() {
  c(
    log_a_L = log(1),
    log_k_C = log(5),
    log_xi_C = log(100),
    log_rho_C = log(1)
  )
}

unpack_pars <- function(x, fixed_pars) {
  c(
    fixed_pars,
    a_L = exp(unname(x["log_a_L"])),
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
    k_C = exp(runif(1, log(1e-4), log(1))),
    xi_C = exp(runif(1, log(1e-3), log(10))),
    rho_C = exp(runif(1, log(1e-4), log(0.5)))
  )
}

make_random_start_set <- function(n_start = 8, seed = 123) {

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

fit_one_parameter_setting_multistart <- function(parameter_row, n_start = 8, seed = 123) {

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
# 12. Fit selected fixed-parameter settings
# ============================================================

selected_fit_list <- vector("list", nrow(fixed_parameter_grid))

for (i in seq_len(nrow(fixed_parameter_grid))) {
  selected_fit_list[[i]] <- fit_one_parameter_setting_multistart(
    parameter_row = fixed_parameter_grid[i, ],
    n_start = N_RANDOM_START,
    seed = RANDOM_SEED_BASE + i
  )
}

fit_list_best <- purrr::map(selected_fit_list, "best_fit")

parameter_fit_table <- purrr::map_dfr(fit_list_best, "summary")

best_fit_params <- purrr::map_dfr(fit_list_best, "params")

all_random_fit_summaries <- purrr::map_dfr(
  seq_along(selected_fit_list),
  function(i) {
    selected_fit_list[[i]]$all_summaries
  }
)

write.csv(
  all_random_fit_summaries,
  file.path(out_dir, "all_random_fit_summaries.csv"),
  row.names = FALSE
)

# ============================================================
# 13. Diagnostics and full outputs for each setting
# ============================================================

setting_diagnostics_list <- vector("list", nrow(fixed_parameter_grid))
v_production_summary_list <- vector("list", nrow(fixed_parameter_grid))
cd8_effect_summary_list <- vector("list", nrow(fixed_parameter_grid))
all_residual_table_list <- vector("list", nrow(fixed_parameter_grid))
all_simulation_list <- vector("list", nrow(fixed_parameter_grid))
all_prediction_list <- vector("list", nrow(fixed_parameter_grid))

for (i in seq_len(nrow(fixed_parameter_grid))) {

  best_fit_i <- fit_list_best[[i]]
  parameter_row_i <- fixed_parameter_grid[i, ]
  fixed_pars_i <- row_to_fixed_pars(parameter_row_i)

  if (is.null(best_fit_i$fit)) {
    next
  }

  residual_table_i <- make_residual_table(
    x = best_fit_i$fit$par,
    fixed_pars = fixed_pars_i
  ) %>%
    mutate(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description
    )

  pars_i <- unpack_pars(best_fit_i$fit$par, fixed_pars_i)

  sim_i <- simulate_all(pars = pars_i) %>%
    mutate(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description
    )

  pred_i <- make_predictions(sim_i) %>%
    mutate(
      setting_id = parameter_row_i$setting_id,
      description = parameter_row_i$description
    )

  all_residual_table_list[[i]] <- residual_table_i
  all_simulation_list[[i]] <- sim_i
  all_prediction_list[[i]] <- pred_i

  setting_diagnostics_list[[i]] <- residual_table_i %>%
    group_by(setting_id, description, variable) %>%
    summarise(
      n = n(),
      n_censored = sum(is_censored, na.rm = TRUE),
      mean_raw_resid = mean(raw_resid, na.rm = TRUE),
      median_raw_resid = median(raw_resid, na.rm = TRUE),
      RMSE_raw = sqrt(mean(raw_resid^2, na.rm = TRUE)),
      weighted_RSS = sum(resid^2, na.rm = TRUE),
      .groups = "drop"
    )

  v_production_summary_list[[i]] <- sim_i %>%
    group_by(setting_id, description, group) %>%
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
    )

  cd8_effect_summary_list[[i]] <- sim_i %>%
    group_by(setting_id, description, id, group) %>%
    summarise(
      C0 = first(C0),
      C0_source = first(C0_source),
      C_base_group = first(C_base_group),
      C_base_source = first(C_base_source),
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
    )
}

all_settings_residual_table <- bind_rows(all_residual_table_list)
all_settings_simulation_full <- bind_rows(all_simulation_list)
all_settings_predictions_long <- bind_rows(all_prediction_list)
setting_diagnostics_by_variable <- bind_rows(setting_diagnostics_list)
v_production_summary_by_setting <- bind_rows(v_production_summary_list)
cd8_effect_summary_by_id <- bind_rows(cd8_effect_summary_list)

write.csv(
  all_settings_residual_table,
  file.path(out_dir, "all_settings_residual_table.csv"),
  row.names = FALSE
)

write.csv(
  all_settings_simulation_full,
  file.path(out_dir, "all_settings_simulation_full.csv"),
  row.names = FALSE
)

write.csv(
  all_settings_predictions_long,
  file.path(out_dir, "all_settings_predictions_long.csv"),
  row.names = FALSE
)

write.csv(
  setting_diagnostics_by_variable,
  file.path(out_dir, "selected_settings_diagnostics_by_variable.csv"),
  row.names = FALSE
)

write.csv(
  v_production_summary_by_setting,
  file.path(out_dir, "v_production_effect_summary_by_setting.csv"),
  row.names = FALSE
)

write.csv(
  cd8_effect_summary_by_id,
  file.path(out_dir, "CD8_effect_summary_by_setting_and_id.csv"),
  row.names = FALSE
)

# ============================================================
# 14. Ranking table
# ============================================================

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

fixed_p_ranking <- parameter_fit_table %>%
  left_join(fixed_parameter_grid, by = c("setting_id", "description")) %>%
  left_join(
    best_fit_params %>% select(setting_id, a_L, p, k_C, xi_C, rho_C),
    by = "setting_id"
  ) %>%
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
    rank_by_problem_score = rank(problem_score, ties.method = "first", na.last = "keep"),
    c_low_flag = c < 10,
    p_low_flag = p < 500,
    xi_C_effective_flag = median_p_eff_over_p_G < 0.95,
    xi_C_strong_suppression_flag = median_p_eff_over_p_G < 0.5,
    practical_note = case_when(
      c_low_flag & p_low_flag ~ "c and fitted p are both low; treat as sensitivity setting",
      c_low_flag ~ "c is low; treat as sensitivity setting",
      p_low_flag ~ "fitted p is low; check biological plausibility",
      TRUE ~ "more suitable as mainline candidate"
    )
  ) %>%
  arrange(problem_score, RSS)

print(fixed_p_ranking)

write.csv(
  fixed_p_ranking,
  file.path(out_dir, "fixed_p_ranking.csv"),
  row.names = FALSE
)

write.csv(
  best_fit_params,
  file.path(out_dir, "fixed_p_best_fit_free_params.csv"),
  row.names = FALSE
)

write.csv(
  parameter_fit_table,
  file.path(out_dir, "fixed_p_fit_summary.csv"),
  row.names = FALSE
)

best_by_problem_score <- fixed_p_ranking %>%
  arrange(problem_score, RSS) %>%
  slice(1)

best_by_RSS <- fixed_p_ranking %>%
  arrange(RSS) %>%
  slice(1)

best_mainline_candidate <- fixed_p_ranking %>%
  filter(!c_low_flag, !p_low_flag) %>%
  arrange(problem_score, RSS) %>%
  slice(1)

if (nrow(best_mainline_candidate) == 0) {
  best_mainline_candidate <- fixed_p_ranking %>%
    arrange(problem_score, RSS) %>%
    slice(1)
}

summary_recommendation <- tibble(
  category = c("best_by_problem_score", "best_by_RSS", "best_mainline_candidate"),
  setting_id = c(
    best_by_problem_score$setting_id,
    best_by_RSS$setting_id,
    best_mainline_candidate$setting_id
  ),
  description = c(
    best_by_problem_score$description,
    best_by_RSS$description,
    best_mainline_candidate$description
  ),
  RSS = c(
    best_by_problem_score$RSS,
    best_by_RSS$RSS,
    best_mainline_candidate$RSS
  ),
  problem_score = c(
    best_by_problem_score$problem_score,
    best_by_RSS$problem_score,
    best_mainline_candidate$problem_score
  ),
  c = c(
    best_by_problem_score$c,
    best_by_RSS$c,
    best_mainline_candidate$c
  ),
  p = c(
    best_by_problem_score$p,
    best_by_RSS$p,
    best_mainline_candidate$p
  ),
  xi_C = c(
    best_by_problem_score$xi_C,
    best_by_RSS$xi_C,
    best_mainline_candidate$xi_C
  ),
  median_p_eff_over_p_G = c(
    best_by_problem_score$median_p_eff_over_p_G,
    best_by_RSS$median_p_eff_over_p_G,
    best_mainline_candidate$median_p_eff_over_p_G
  ),
  practical_note = c(
    best_by_problem_score$practical_note,
    best_by_RSS$practical_note,
    best_mainline_candidate$practical_note
  )
)

print(summary_recommendation)

write.csv(
  summary_recommendation,
  file.path(out_dir, "summary_recommendation.csv"),
  row.names = FALSE
)

# ============================================================
# 15. Save best-mainline detailed files
# ============================================================

best_mainline_setting_id <- best_mainline_candidate$setting_id

best_mainline_residual_table <- all_settings_residual_table %>%
  filter(setting_id == best_mainline_setting_id)

best_mainline_simulation_full <- all_settings_simulation_full %>%
  filter(setting_id == best_mainline_setting_id)

best_mainline_predictions_long <- all_settings_predictions_long %>%
  filter(setting_id == best_mainline_setting_id)

write.csv(
  best_mainline_residual_table,
  file.path(out_dir, "best_mainline_residual_table.csv"),
  row.names = FALSE
)

write.csv(
  best_mainline_simulation_full,
  file.path(out_dir, "best_mainline_simulation_full.csv"),
  row.names = FALSE
)

write.csv(
  best_mainline_predictions_long,
  file.path(out_dir, "best_mainline_predictions_long.csv"),
  row.names = FALSE
)

# ============================================================
# 16. Figures
# ============================================================

p_rank_problem <- fixed_p_ranking %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  mutate(setting_label = factor(setting_label, levels = rev(setting_label))) %>%
  ggplot(aes(x = setting_label, y = problem_score)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Fixed p values ranked by problem score",
    x = "Fixed-parameter setting",
    y = "Problem score = weighted RSS + pVL/CD8 bias penalty"
  )

print(p_rank_problem)
save_plot_png_pdf(p_rank_problem, "Fig1_fixed_p_ranking_problem_score", width = 9, height = 5)

p_rank_RSS <- fixed_p_ranking %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  mutate(setting_label = factor(setting_label, levels = rev(setting_label))) %>%
  ggplot(aes(x = setting_label, y = RSS)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Fixed p values ranked by weighted RSS",
    x = "Fixed-parameter setting",
    y = "Weighted RSS"
  )

print(p_rank_RSS)
save_plot_png_pdf(p_rank_RSS, "Fig2_fixed_p_ranking_RSS", width = 9, height = 5)

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
    title = "Mean raw residual bias for pVL and CD8",
    x = "Fixed-parameter setting",
    y = "Mean raw residual: predicted - observed"
  )

print(p_bias_compare)
save_plot_png_pdf(p_bias_compare, "Fig3_pVL_CD8_bias_selected_settings", width = 10, height = 5)

p_xiC_compare <- v_production_summary_by_setting %>%
  mutate(setting_label = paste(setting_id, description, sep = ": ")) %>%
  ggplot(aes(x = setting_label, y = median_p_eff_over_p_G)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ group) +
  theme_bw() +
  labs(
    title = "Effective viral production fraction under xi_C",
    x = "Fixed-parameter setting",
    y = "Median p_eff / p"
  )

print(p_xiC_compare)
save_plot_png_pdf(p_xiC_compare, "Fig4_p_eff_over_p_selected_settings", width = 10, height = 5)

p_obs_pred_all <- all_settings_residual_table %>%
  ggplot(aes(x = value, y = pred)) +
  geom_point(aes(shape = is_censored), alpha = 0.65, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(variable ~ setting_id, scales = "free") +
  theme_bw() +
  labs(
    title = "Observed vs predicted across selected fixed-parameter settings",
    x = "Observed value used for fitting",
    y = "Predicted value",
    shape = "Censored"
  )

print(p_obs_pred_all)
save_plot_png_pdf(p_obs_pred_all, "Fig5_observed_vs_predicted_all_selected_settings", width = 12, height = 8)

p_timecourse_best <- best_mainline_residual_table %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = pred, group = id), alpha = 0.35, linewidth = 0.5) +
  geom_point(aes(y = value, shape = is_censored), alpha = 0.7, size = 1.8) +
  facet_grid(variable ~ group, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste0("Observed and predicted time courses: mainline candidate ", best_mainline_setting_id),
    x = "Time",
    y = "Value used for fitting",
    shape = "Censored"
  )

print(p_timecourse_best)
save_plot_png_pdf(p_timecourse_best, "Fig6_timecourse_best_mainline_candidate", width = 10, height = 8)

p_cd8_best <- best_mainline_simulation_full %>%
  select(
    id,
    group,
    time,
    C,
    C_scaled,
    C_base_group,
    CD8_activation_fraction,
    CD8_killing_saturation,
    CD8_killing_rate,
    p_eff_over_p_G
  ) %>%
  pivot_longer(
    cols = c(C, C_scaled, C_base_group, CD8_activation_fraction, CD8_killing_saturation, CD8_killing_rate, p_eff_over_p_G),
    names_to = "quantity",
    values_to = "value"
  ) %>%
  ggplot(aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  facet_grid(quantity ~ group, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste0("CD8 and viral-production effects: mainline candidate ", best_mainline_setting_id),
    x = "Time",
    y = "Value"
  )

print(p_cd8_best)
save_plot_png_pdf(p_cd8_best, "Fig7_CD8_and_xiC_effects_best_mainline_candidate", width = 11, height = 9)

# ============================================================
# 17. Save run settings and final message
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
    "tested_fixed_p_values_under_S03_background",
    "CD8_equation",
    "xi_C_structure",
    "fitted_common_parameters",
    "ranking_rule",
    "mainline_candidate_rule",
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
    paste(fixed_parameter_grid$setting_id, collapse = ","),
    "dC = d_C * (C_base_group - C) + rho_C * I/(K_C + I) * C",
    "p_eff = p / (1 + xi_C * C_scaled), C_scaled = C / 1e6",
    "a_L, k_C, xi_C, rho_C; p is fixed",
    "problem_score = RSS + penalty for pVL underprediction and CD8 overprediction",
    "Among settings without c_low_flag and p_low_flag, choose lowest problem_score; if none, choose lowest problem_score overall",
    out_dir
  )
)

write.csv(
  run_settings,
  file.path(out_dir, "run_settings.csv"),
  row.names = FALSE
)

cat("\nFixed-p scan finished.\n")

cat("\nBest by weighted RSS:\n")
print(best_by_RSS %>% select(setting_id, description, RSS, problem_score, c, p, xi_C, practical_note))

cat("\nBest by problem score:\n")
print(best_by_problem_score %>% select(setting_id, description, RSS, problem_score, c, p, xi_C, practical_note))

cat("\nBest mainline candidate considering practical flags:\n")
print(best_mainline_candidate %>% select(setting_id, description, RSS, problem_score, c, p, xi_C, median_p_eff_over_p_G, practical_note))

cat("\nMain files to check:\n")
cat("1. fixed_p_ranking.csv\n")
cat("2. summary_recommendation.csv\n")
cat("3. selected_settings_diagnostics_by_variable.csv\n")
cat("4. v_production_effect_summary_by_setting.csv\n")
cat("5. CD8_effect_summary_by_setting_and_id.csv\n")
cat("6. all_settings_residual_table.csv\n")
cat("7. all_settings_simulation_full.csv\n")
cat("8. best_mainline_residual_table.csv\n")
cat("9. best_mainline_simulation_full.csv\n")

cat("\nMain figures to check:\n")
cat("1. Fig1_fixed_p_ranking_problem_score.png\n")
cat("2. Fig2_fixed_p_ranking_RSS.png\n")
cat("3. Fig3_pVL_CD8_bias_selected_settings.png\n")
cat("4. Fig4_p_eff_over_p_selected_settings.png\n")
cat("5. Fig5_observed_vs_predicted_all_selected_settings.png\n")
cat("6. Fig6_timecourse_best_mainline_candidate.png\n")
cat("7. Fig7_CD8_and_xiC_effects_best_mainline_candidate.png\n")

cat("\nAll results saved to:\n")
cat(out_dir, "\n")

cat("\nFigures saved to:\n")
cat(fig_dir, "\n")
