library(deSolve)
library(tidyverse)
library(minpack.lm)

# ============================================================
# Total-CD8 saturated-killing model selection script
# Fixed-p version; four-mechanism early/late model selection
# ============================================================
# Current model:
#   T, L, I, V, C
#   C = total CD8 state variable
#
# Fixed parameters based on the previous S03 + fixed-p sensitivity analysis:
#   p       = 8000
#   delta0  = 0.45
#   c       = 23
#   d_C     = 0.75
#   K_C     = 1e4
#   K_kill  = 5
#
# Parameters fitted as common baseline values:
#   a_L, k_C, xi_C, rho_C
#
# Model selection tests whether early vs late groups differ in:
#   1) a_L    : latent-cell activation rate
#   2) k_C    : maximal saturated CD8 cytolytic killing rate
#   3) xi_C   : CD8 non-cytolytic suppression coefficient
#   4) rho_C  : CD8 proliferation / activation response
#
# Sign convention:
#   parameter_early = parameter_late * exp(theta)
#   theta > 0 means early > late.
#
# Important structural update:
#   Old CD8 equation:
#     dC = d_C * (C0 - C) + rho_C_G * I/(K_C + I) * C
#
#   New CD8 equation:
#     dC = d_C * (C_base_group - C) + rho_C_G * I/(K_C + I) * C
#
#   C0 is now only the initial CD8 state.
#   C_base_group is the group-level CD8 homeostatic baseline.
# ============================================================

# ============================================================
# 0. Global settings
# ============================================================

EPS <- 1e-12

# Stage 1 uses broad random exploration.
# Stage 2 reuses Stage 1 best solutions as anchors for all models.
# If runtime is too long, reduce these two numbers first.
N_RANDOM_START_STAGE1 <- 10
N_RANDOM_START_STAGE2 <- 5
N_RANDOM_START <- N_RANDOM_START_STAGE1
RUN_TWO_STAGE_ANCHORED_FITTING <- TRUE
MAX_STAGE2_ANCHORS_PER_MODEL <- 32

RANDOM_SEED_BASE <- 1000

USE_VARIABLE_SCALING <- TRUE

FIT_VARIABLES <- c("pVL", "DNA", "CD4", "CD8")

RUN_SINGLE_SIM_TEST <- TRUE
RUN_DEBUG_M0 <- TRUE
RUN_TEST_M0 <- TRUE

# Fixed constants from S03 background + fixed-p sensitivity analysis.
P_FIXED <- 8000
DELTA0_FIXED <- 0.45
C_FIXED <- 23
D_C_FIXED <- 0.75
K_C_FIXED <- 1e4
K_KILL_FIXED <- 5.0
CD8_SCALE <- 1e6

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ============================================================
# Path settings
# ============================================================

DATA_DIR <- "D:/R-4.6.0/R-ATI-project/Rdata"

DATA_FILE <- file.path(DATA_DIR, "raw_data_with_cd8_counts.csv")
INIT_FILE <- file.path(DATA_DIR, "init_table_total_CD8_no_exhaustion.csv")

out_dir <- file.path(
  DATA_DIR,
  "model_selection_fixed_p8000_S03_CD8_Cbase_four_mechanisms",
  run_tag
)

fig_dir <- file.path(out_dir, "figures")

cat("\nDATA_DIR:\n")
cat(DATA_DIR, "\n")

cat("\nDATA_FILE:\n")
cat(DATA_FILE, "\n")
cat("DATA_FILE exists: ", file.exists(DATA_FILE), "\n")

cat("\nINIT_FILE:\n")
cat(INIT_FILE, "\n")
cat("INIT_FILE exists: ", file.exists(INIT_FILE), "\n")

if (!file.exists(DATA_FILE)) {
  stop("找不到 DATA_FILE: ", DATA_FILE, "\n当前工作目录是: ", getwd())
}

if (!file.exists(INIT_FILE)) {
  stop("找不到 INIT_FILE: ", INIT_FILE, "\n当前工作目录是: ", getwd())
}

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 1. ODE model
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
    xi_C_G <- max(xi_C_G, EPS)

    C_scaled <- C / CD8_scale

    CD8_activation_fraction <- I / (K_C + I)
    CD8_killing_saturation <- C_scaled / (K_kill + C_scaled)
    CD8_killing_rate <- k_C_G * CD8_killing_saturation

    p_eff <- p / (1 + xi_C_G * C_scaled)

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

    dC <- d_C * (C_base_group - C) +
      rho_C_G * CD8_activation_fraction * C

    list(c(dT, dL, dI, dV, dC))
  })
}

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
# 3. Fixed parameters
# ============================================================

fixed_pars <- c(
  d_T = 0.01,
  beta = 1.58e-8,
  f_L = 0.005,
  d_L = 5e-4,
  delta0 = DELTA0_FIXED,
  c = C_FIXED,
  d_C = D_C_FIXED,
  K_C = K_C_FIXED,
  K_kill = K_KILL_FIXED,
  CD8_scale = CD8_SCALE,
  p = P_FIXED
)

# ============================================================
# 4. Default free-parameter starts
# ============================================================
# p is fixed and is NOT fitted.
# a_L, k_C, xi_C, rho_C are fitted as baseline late-group values.
# theta_* controls early/late ratios when the corresponding model flag is on.

default_start_free_pars <- c(
  a_L = 1e-3,
  k_C = 0.01,
  xi_C = 1,
  rho_C = 0.05,
  theta_a = 0,
  theta_k = 0,
  theta_xi = 0,
  theta_rho = 0
)

# ============================================================
# 5. Test parameters
# ============================================================

test_pars <- c(
  fixed_pars,
  a_L = 1e-3,
  k_C = 0.01,
  xi_C = 1,
  rho_C = 0.05,
  theta_a = 0,
  theta_k = 0,
  theta_xi = 0,
  theta_rho = 0
)

test_flags <- c(
  za = 0,
  zk = 0,
  zxi = 0,
  zrho = 0
)

# ============================================================
# 6. Helper functions
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

safe_log10 <- function(x) {
  log10(pmax(x, EPS))
}

# ============================================================
# 7. Data cleaning
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
  filter(
    is.na(T0) | is.na(L0) | is.na(I0) | is.na(V0)
  )

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
  filter(
    !is.na(T0),
    !is.na(L0),
    !is.na(I0),
    !is.na(V0)
  )

# ============================================================
# 8. Data scale checks
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
# 9. Add initial CD8 and group-level CD8 baseline
# ============================================================
# C0 is only the initial state.
# C_base_group is the CD8 homeostatic set point in the new CD8 equation.
#
# Priority for C0:
#   1) CD8_0 already provided in init table
#   2) individual CD8 at time 0
#   3) individual median CD8
#   4) group median CD8
#   5) overall median CD8
#
# C_base_group:
#   group-level median CD8 across all available CD8 observations.

cd8_data <- data_raw %>%
  filter(
    variable == "CD8",
    !is.na(value),
    is.finite(value)
  )

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
  left_join(
    cd8_by_id,
    by = c("id", "group")
  ) %>%
  left_join(
    cd8_by_group,
    by = "group"
  ) %>%
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
# 10. Fitting data preprocessing
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
  filter(
    !is.na(value),
    is.finite(value)
  )

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
# 11. Residual scaling
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
# 12. Candidate models
# ============================================================
# za   : a_L early/late difference
# zk   : k_C early/late difference
# zxi  : xi_C early/late difference
# zrho : rho_C early/late difference
#
# p is fixed in all candidate models.
# Exhaustion compartment is removed.

make_model_name <- function(za, zk, zxi, zrho) {

  labels <- c(
    ifelse(za == 1, "a", NA_character_),
    ifelse(zk == 1, "k", NA_character_),
    ifelse(zxi == 1, "xi", NA_character_),
    ifelse(zrho == 1, "rho", NA_character_)
  )

  labels <- labels[!is.na(labels)]

  if (length(labels) == 0) {
    return("M0")
  }

  paste0("M", paste(labels, collapse = "_"))
}

candidate_models <- expand.grid(
  za = c(0, 1),
  zk = c(0, 1),
  zxi = c(0, 1),
  zrho = c(0, 1)
) %>%
  as_tibble() %>%
  mutate(
    n_group_effects = za + zk + zxi + zrho,
    model = pmap_chr(
      list(za, zk, zxi, zrho),
      make_model_name
    )
  ) %>%
  arrange(n_group_effects, za, zk, zxi, zrho) %>%
  select(model, za, zk, zxi, zrho)

print(candidate_models)

write.csv(
  candidate_models,
  file.path(out_dir, "candidate_models.csv"),
  row.names = FALSE
)

# ============================================================
# 13. Single-subject simulation
# ============================================================

simulate_one <- function(id_row, pars, times, model_flags) {

  G <- ifelse(is_early_group(id_row$group), 1, 0)

  a_L_G <- as.numeric(pars["a_L"]) *
    exp(as.numeric(model_flags["za"]) * as.numeric(pars["theta_a"]) * G)

  k_C_G <- as.numeric(pars["k_C"]) *
    exp(as.numeric(model_flags["zk"]) * as.numeric(pars["theta_k"]) * G)

  xi_C_G <- as.numeric(pars["xi_C"]) *
    exp(as.numeric(model_flags["zxi"]) * as.numeric(pars["theta_xi"]) * G)

  rho_C_G <- as.numeric(pars["rho_C"]) *
    exp(as.numeric(model_flags["zrho"]) * as.numeric(pars["theta_rho"]) * G)

  T0 <- max(as.numeric(id_row$T0), EPS)
  L0 <- max(as.numeric(id_row$L0), EPS)
  I0 <- max(as.numeric(id_row$I0), EPS)
  V0 <- max(as.numeric(id_row$V0), EPS)
  C0 <- max(as.numeric(id_row$C0), EPS)
  C_base_group <- max(as.numeric(id_row$C_base_group), EPS)

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
        a_L_G = a_L_G,
        k_C_G = k_C_G,
        xi_C_G = xi_C_G,
        rho_C_G = rho_C_G
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
      G_early = G,
      T0 = T0,
      L0 = L0,
      I0 = I0,
      V0 = V0,
      C0 = C0,
      C0_source = id_row$C0_source,
      C_base_group = C_base_group,
      C_base_source = id_row$C_base_source,
      a_L_G = a_L_G,
      p_G = as.numeric(pars["p"]),
      k_C_G = k_C_G,
      xi_C_G = xi_C_G,
      rho_C_G = rho_C_G,
      C_scaled = C / as.numeric(pars["CD8_scale"]),
      CD8_activation_fraction = I / (as.numeric(pars["K_C"]) + I),
      CD8_killing_saturation = C_scaled / (as.numeric(pars["K_kill"]) + C_scaled),
      CD8_killing_rate = k_C_G * CD8_killing_saturation,
      p_eff = as.numeric(pars["p"]) / (1 + xi_C_G * C_scaled),
      p_eff_over_p_G = p_eff / as.numeric(pars["p"]),
      CD4_total = T + L + I,
      CD8_total = C,
      pred_pVL_log10 = safe_log10(V),
      pred_DNA_log10 = safe_log10(1e6 * (L + I) / pmax(CD4_total, EPS)),
      pred_CD4_log10 = safe_log10(CD4_total),
      pred_CD8_log10 = safe_log10(CD8_total)
    )
}

# ============================================================
# 14. Free-parameter vector
# ============================================================

make_start_vector <- function(model_flags, start_free_pars = default_start_free_pars) {

  x0 <- c(
    log_a_L = log(unname(start_free_pars["a_L"])),
    log_k_C = log(unname(start_free_pars["k_C"])),
    log_xi_C = log(unname(start_free_pars["xi_C"])),
    log_rho_C = log(unname(start_free_pars["rho_C"]))
  )

  if (model_flags["za"] == 1) {
    x0 <- c(x0, theta_a = unname(start_free_pars["theta_a"]))
  }

  if (model_flags["zk"] == 1) {
    x0 <- c(x0, theta_k = unname(start_free_pars["theta_k"]))
  }

  if (model_flags["zxi"] == 1) {
    x0 <- c(x0, theta_xi = unname(start_free_pars["theta_xi"]))
  }

  if (model_flags["zrho"] == 1) {
    x0 <- c(x0, theta_rho = unname(start_free_pars["theta_rho"]))
  }

  x0
}

make_lower_vector <- function(model_flags) {

  lower <- c(
    log_a_L = log(1e-8),
    log_k_C = log(1e-6),
    log_xi_C = log(1e-5),
    log_rho_C = log(1e-6)
  )

  if (model_flags["za"] == 1) {
    lower <- c(lower, theta_a = -5)
  }

  if (model_flags["zk"] == 1) {
    lower <- c(lower, theta_k = -5)
  }

  if (model_flags["zxi"] == 1) {
    lower <- c(lower, theta_xi = -5)
  }

  if (model_flags["zrho"] == 1) {
    lower <- c(lower, theta_rho = -5)
  }

  lower
}

make_upper_vector <- function(model_flags) {

  upper <- c(
    log_a_L = log(1),
    log_k_C = log(10),
    log_xi_C = log(100),
    log_rho_C = log(5)
  )

  if (model_flags["za"] == 1) {
    upper <- c(upper, theta_a = 5)
  }

  if (model_flags["zk"] == 1) {
    upper <- c(upper, theta_k = 5)
  }

  if (model_flags["zxi"] == 1) {
    upper <- c(upper, theta_xi = 5)
  }

  if (model_flags["zrho"] == 1) {
    upper <- c(upper, theta_rho = 5)
  }

  upper
}

unpack_pars <- function(x, model_flags) {

  pars <- c(
    fixed_pars,
    a_L = exp(unname(x["log_a_L"])),
    k_C = exp(unname(x["log_k_C"])),
    xi_C = exp(unname(x["log_xi_C"])),
    rho_C = exp(unname(x["log_rho_C"])),
    theta_a = 0,
    theta_k = 0,
    theta_xi = 0,
    theta_rho = 0
  )

  if (model_flags["za"] == 1) {
    pars["theta_a"] <- unname(x["theta_a"])
  }

  if (model_flags["zk"] == 1) {
    pars["theta_k"] <- unname(x["theta_k"])
  }

  if (model_flags["zxi"] == 1) {
    pars["theta_xi"] <- unname(x["theta_xi"])
  }

  if (model_flags["zrho"] == 1) {
    pars["theta_rho"] <- unname(x["theta_rho"])
  }

  pars
}

# ============================================================
# 15. Simulate all individuals
# ============================================================

simulate_all <- function(pars, model_flags) {

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
      times = times_i,
      model_flags = model_flags
    )

    if (is.null(sim_i)) {
      return(NULL)
    }

    sim_list[[i]] <- sim_i
  }

  bind_rows(sim_list)
}

# ============================================================
# 16. Predictions
# ============================================================

make_predictions <- function(sim) {

  pred_cols <- c(
    pred_pVL_log10 = "pVL",
    pred_DNA_log10 = "DNA",
    pred_CD4_log10 = "CD4",
    pred_CD8_log10 = "CD8"
  )

  selected_pred_cols <- names(pred_cols)[pred_cols %in% FIT_VARIABLES]

  sim %>%
    select(
      id,
      group,
      time,
      all_of(selected_pred_cols)
    ) %>%
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
# 17. Residuals
# ============================================================

make_residual_table <- function(x, model_flags) {

  pars <- unpack_pars(x, model_flags)

  sim <- simulate_all(
    pars = pars,
    model_flags = model_flags
  )

  if (is.null(sim)) {
    return(NULL)
  }

  pred <- make_predictions(sim)

  joined <- data_fit %>%
    left_join(
      pred,
      by = c("id", "group", "time", "variable")
    ) %>%
    left_join(
      residual_scale_table,
      by = "variable"
    )

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

make_residuals <- function(x, model_flags) {

  joined <- make_residual_table(
    x = x,
    model_flags = model_flags
  )

  if (is.null(joined)) {
    return(rep(1e6, nrow(data_fit)))
  }

  resid <- joined$resid

  if (any(!is.finite(resid))) {
    return(rep(1e6, nrow(data_fit)))
  }

  resid
}

# ============================================================
# 18. Debug residuals
# ============================================================

debug_residuals <- function(x, model_flags) {

  pars <- unpack_pars(x, model_flags)

  cat("\n--- Parameter check ---\n")
  print(pars)

  sim <- simulate_all(
    pars = pars,
    model_flags = model_flags
  )

  if (is.null(sim)) {
    cat("\nSimulation returned NULL.\n")
    return(invisible(NULL))
  }

  cat("\nSimulation succeeded.\n")
  cat("sim rows:", nrow(sim), "\n")

  cat("\nCheck simulated values:\n")
  print(
    sim %>%
      summarise(
        min_T = min(T, na.rm = TRUE),
        max_T = max(T, na.rm = TRUE),
        min_L = min(L, na.rm = TRUE),
        max_L = max(L, na.rm = TRUE),
        min_I = min(I, na.rm = TRUE),
        max_I = max(I, na.rm = TRUE),
        min_V = min(V, na.rm = TRUE),
        max_V = max(V, na.rm = TRUE),
        min_C = min(C, na.rm = TRUE),
        max_C = max(C, na.rm = TRUE),
        min_CD8_killing_rate = min(CD8_killing_rate, na.rm = TRUE),
        max_CD8_killing_rate = max(CD8_killing_rate, na.rm = TRUE),
        min_p_eff = min(p_eff, na.rm = TRUE),
        max_p_eff = max(p_eff, na.rm = TRUE),
        any_nonfinite = any(
          !is.finite(T) |
            !is.finite(L) |
            !is.finite(I) |
            !is.finite(V) |
            !is.finite(C)
        )
      )
  )

  pred <- make_predictions(sim)

  joined <- data_fit %>%
    left_join(
      pred,
      by = c("id", "group", "time", "variable")
    ) %>%
    left_join(
      residual_scale_table,
      by = "variable"
    )

  cat("\nJoined rows:", nrow(joined), "\n")
  cat("Number of missing predictions:", sum(is.na(joined$pred)), "\n")
  cat("Number of non-finite predictions:", sum(!is.finite(joined$pred)), "\n")

  if (any(is.na(joined$pred))) {
    cat("\nRows with missing predictions:\n")
    print(
      joined %>%
        filter(is.na(pred)) %>%
        select(id, group, time, variable, value) %>%
        head(30)
    )
  }

  joined <- joined %>%
    mutate(
      raw_resid = case_when(
        is_censored & pred <= lod ~ 0,
        is_censored & pred > lod ~ pred - lod,
        TRUE ~ pred - value
      ),
      resid = raw_resid / residual_scale
    )

  cat("\nResidual check:\n")
  cat("Weighted RSS:", sum(joined$resid^2, na.rm = TRUE), "\n")
  cat("Raw RSS:", sum(joined$raw_resid^2, na.rm = TRUE), "\n")
  cat("Any non-finite residual:", any(!is.finite(joined$resid)), "\n")

  cat("\nResidual summary by variable:\n")
  print(
    joined %>%
      group_by(variable) %>%
      summarise(
        n = n(),
        n_censored = sum(is_censored, na.rm = TRUE),
        raw_RSS = sum(raw_resid^2, na.rm = TRUE),
        weighted_RSS = sum(resid^2, na.rm = TRUE),
        mean_abs_raw_resid = mean(abs(raw_resid), na.rm = TRUE),
        .groups = "drop"
      )
  )

  invisible(joined)
}

# ============================================================
# 19. Fit one candidate model
# ============================================================

fit_one_model <- function(
    model_row,
    start_free_pars = default_start_free_pars,
    start_id = NA_integer_
) {

  model_flags <- c(
    za = as.numeric(model_row$za),
    zk = as.numeric(model_row$zk),
    zxi = as.numeric(model_row$zxi),
    zrho = as.numeric(model_row$zrho)
  )

  x0 <- make_start_vector(
    model_flags = model_flags,
    start_free_pars = start_free_pars
  )

  fit <- tryCatch({
    nls.lm(
      par = x0,
      lower = make_lower_vector(model_flags),
      upper = make_upper_vector(model_flags),
      fn = make_residuals,
      model_flags = model_flags,
      control = nls.lm.control(
        maxiter = 300,
        ftol = 1e-10,
        ptol = 1e-10,
        gtol = 1e-10
      )
    )
  }, error = function(e) {
    message("Fitting failed for model ", model_row$model, ": ", e$message)
    return(NULL)
  })

  failed_summary <- function(convergence_info = NA_integer_) {
    tibble(
      model = model_row$model,
      za = model_row$za,
      zk = model_row$zk,
      zxi = model_row$zxi,
      zrho = model_row$zrho,
      n = nrow(data_fit),
      q = length(x0),
      K = length(x0) + 1,
      RSS = NA_real_,
      raw_RSS = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      convergence_info = convergence_info,
      start_id = start_id
    )
  }

  failed_params <- function() {
    tibble(
      model = model_row$model,
      a_L = NA_real_,
      k_C = NA_real_,
      xi_C = NA_real_,
      rho_C = NA_real_,
      theta_a = NA_real_,
      theta_k = NA_real_,
      theta_xi = NA_real_,
      theta_rho = NA_real_,
      a_L_late = NA_real_,
      a_L_early = NA_real_,
      k_C_late = NA_real_,
      k_C_early = NA_real_,
      xi_C_late = NA_real_,
      xi_C_early = NA_real_,
      rho_C_late = NA_real_,
      rho_C_early = NA_real_,
      ratio_a_L_early_late = NA_real_,
      ratio_k_C_early_late = NA_real_,
      ratio_xi_C_early_late = NA_real_,
      ratio_rho_C_early_late = NA_real_,
      fixed_p = P_FIXED,
      start_id = start_id
    )
  }

  if (is.null(fit)) {
    return(list(
      fit = NULL,
      summary = failed_summary(),
      params = failed_params()
    ))
  }

  x_hat <- fit$par
  pars_hat <- unpack_pars(x_hat, model_flags)

  residual_table_hat <- make_residual_table(
    x = x_hat,
    model_flags = model_flags
  )

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
  RSS_safe <- max(RSS, .Machine$double.eps)

  n <- length(resid_hat)
  q <- length(x_hat)
  K <- q + 1

  logLik <- -0.5 * n * (log(2 * pi) + 1 + log(RSS_safe / n))

  AIC <- -2 * logLik + 2 * K

  AICc <- ifelse(
    n > K + 1,
    AIC + (2 * K * (K + 1)) / (n - K - 1),
    NA_real_
  )

  a_L_late <- unname(pars_hat["a_L"])
  a_L_early <- unname(pars_hat["a_L"]) *
    exp(unname(model_flags["za"]) * unname(pars_hat["theta_a"]))

  k_C_late <- unname(pars_hat["k_C"])
  k_C_early <- unname(pars_hat["k_C"]) *
    exp(unname(model_flags["zk"]) * unname(pars_hat["theta_k"]))

  xi_C_late <- unname(pars_hat["xi_C"])
  xi_C_early <- unname(pars_hat["xi_C"]) *
    exp(unname(model_flags["zxi"]) * unname(pars_hat["theta_xi"]))

  rho_C_late <- unname(pars_hat["rho_C"])
  rho_C_early <- unname(pars_hat["rho_C"]) *
    exp(unname(model_flags["zrho"]) * unname(pars_hat["theta_rho"]))

  summary_row <- tibble(
    model = model_row$model,
    za = model_row$za,
    zk = model_row$zk,
    zxi = model_row$zxi,
    zrho = model_row$zrho,
    n = n,
    q = q,
    K = K,
    RSS = RSS,
    raw_RSS = raw_RSS,
    logLik = logLik,
    AIC = AIC,
    AICc = AICc,
    convergence_info = fit$info,
    start_id = start_id
  )

  param_row <- tibble(
    model = model_row$model,
    a_L = unname(pars_hat["a_L"]),
    k_C = unname(pars_hat["k_C"]),
    xi_C = unname(pars_hat["xi_C"]),
    rho_C = unname(pars_hat["rho_C"]),
    theta_a = unname(pars_hat["theta_a"]),
    theta_k = unname(pars_hat["theta_k"]),
    theta_xi = unname(pars_hat["theta_xi"]),
    theta_rho = unname(pars_hat["theta_rho"]),
    a_L_late = a_L_late,
    a_L_early = a_L_early,
    k_C_late = k_C_late,
    k_C_early = k_C_early,
    xi_C_late = xi_C_late,
    xi_C_early = xi_C_early,
    rho_C_late = rho_C_late,
    rho_C_early = rho_C_early,
    ratio_a_L_early_late = a_L_early / a_L_late,
    ratio_k_C_early_late = k_C_early / k_C_late,
    ratio_xi_C_early_late = xi_C_early / xi_C_late,
    ratio_rho_C_early_late = rho_C_early / rho_C_late,
    fixed_p = P_FIXED,
    start_id = start_id
  )

  list(
    fit = fit,
    summary = summary_row,
    params = param_row
  )
}

# ============================================================
# 20. Random starts
# ============================================================

make_random_start <- function() {

  c(
    a_L = exp(runif(1, log(1e-5), log(1e-1))),
    k_C = exp(runif(1, log(1e-4), log(1))),
    xi_C = exp(runif(1, log(1e-3), log(10))),
    rho_C = exp(runif(1, log(1e-4), log(1))),
    theta_a = runif(1, -2, 2),
    theta_k = runif(1, -2, 2),
    theta_xi = runif(1, -2, 2),
    theta_rho = runif(1, -2, 2)
  )
}

make_random_start_set <- function(
    n_start = 10,
    seed = 123,
    anchor_starts = list()
) {

  set.seed(seed)

  start_set <- list()

  start_set[[length(start_set) + 1]] <- default_start_free_pars

  if (length(anchor_starts) > 0) {
    for (i in seq_along(anchor_starts)) {
      start_set[[length(start_set) + 1]] <- anchor_starts[[i]]
    }
  }

  for (i in seq_len(n_start)) {
    start_set[[length(start_set) + 1]] <- make_random_start()
  }

  start_set
}

fit_one_model_random_multistart <- function(
    model_row,
    n_start = 10,
    seed = 123,
    anchor_starts = list()
) {

  start_set <- make_random_start_set(
    n_start = n_start,
    seed = seed,
    anchor_starts = anchor_starts
  )

  fit_candidates <- purrr::map(
    seq_along(start_set),
    function(i) {
      cat(
        "\nModel:", model_row$model,
        "| start:", i,
        "\n"
      )

      fit_one_model(
        model_row = model_row,
        start_free_pars = start_set[[i]],
        start_id = i
      )
    }
  )

  summary_candidates <- purrr::map_dfr(
    fit_candidates,
    "summary"
  ) %>%
    mutate(
      n_random_start = n_start,
      n_anchor_start = length(anchor_starts)
    )

  valid_candidates <- summary_candidates %>%
    filter(!is.na(logLik), is.finite(logLik))

  if (nrow(valid_candidates) == 0) {

    failed_summary <- summary_candidates %>%
      slice(1) %>%
      mutate(
        best_start_id = NA_integer_
      )

    failed_params <- purrr::map_dfr(
      fit_candidates,
      "params"
    ) %>%
      slice(1) %>%
      mutate(
        best_start_id = NA_integer_,
        n_random_start = n_start,
        n_anchor_start = length(anchor_starts)
      )

    return(list(
      best_fit = list(
        fit = NULL,
        summary = failed_summary,
        params = failed_params
      ),
      all_summaries = summary_candidates
    ))
  }

  best_start_id <- valid_candidates %>%
    arrange(desc(logLik)) %>%
    slice(1) %>%
    pull(start_id)

  best_fit <- fit_candidates[[best_start_id]]

  best_fit$summary <- best_fit$summary %>%
    mutate(
      best_start_id = best_start_id,
      n_random_start = n_start,
      n_anchor_start = length(anchor_starts)
    )

  best_fit$params <- best_fit$params %>%
    mutate(
      best_start_id = best_start_id,
      n_random_start = n_start,
      n_anchor_start = length(anchor_starts)
    )

  list(
    best_fit = best_fit,
    all_summaries = summary_candidates
  )
}

# ============================================================
# 21. Anchor starts from nested parent models
# ============================================================

params_row_to_start_free <- function(param_row) {

  c(
    a_L = as.numeric(param_row$a_L),
    k_C = as.numeric(param_row$k_C),
    xi_C = as.numeric(param_row$xi_C),
    rho_C = as.numeric(param_row$rho_C),
    theta_a = ifelse(is.na(param_row$theta_a), 0, as.numeric(param_row$theta_a)),
    theta_k = ifelse(is.na(param_row$theta_k), 0, as.numeric(param_row$theta_k)),
    theta_xi = ifelse(is.na(param_row$theta_xi), 0, as.numeric(param_row$theta_xi)),
    theta_rho = ifelse(is.na(param_row$theta_rho), 0, as.numeric(param_row$theta_rho))
  )
}

is_nested_parent <- function(parent_row, child_row) {

  parent_flags <- c(
    za = as.numeric(parent_row$za),
    zk = as.numeric(parent_row$zk),
    zxi = as.numeric(parent_row$zxi),
    zrho = as.numeric(parent_row$zrho)
  )

  child_flags <- c(
    za = as.numeric(child_row$za),
    zk = as.numeric(child_row$zk),
    zxi = as.numeric(child_row$zxi),
    zrho = as.numeric(child_row$zrho)
  )

  all(parent_flags <= child_flags) &&
    any(parent_flags < child_flags)
}

# ============================================================
# 22. Optional single-subject simulation test
# ============================================================

if (RUN_SINGLE_SIM_TEST) {

  times <- seq(
    0,
    max(data_fit$time, na.rm = TRUE),
    by = 0.1
  )

  test_id_row <- init_fit %>%
    arrange(desc(L0)) %>%
    slice(1)

  sim_test <- simulate_one(
    id_row = test_id_row,
    pars = test_pars,
    times = times,
    model_flags = test_flags
  )

  if (is.null(sim_test)) {
    stop("Simulation failed.")
  }

  cat("\nSingle-subject total-CD8 saturated-killing fixed-p model check:\n")
  print(
    sim_test %>%
      summarise(
        id = first(id),
        group = first(group),
        C0 = first(C0),
        C_base_group = first(C_base_group),
        median_C = median(C, na.rm = TRUE),
        median_CD8_activation_fraction = median(CD8_activation_fraction, na.rm = TRUE),
        median_CD8_killing_saturation = median(CD8_killing_saturation, na.rm = TRUE),
        median_CD8_killing_rate = median(CD8_killing_rate, na.rm = TRUE),
        median_p_eff_over_p_G = median(p_eff_over_p_G, na.rm = TRUE),
        L0 = first(L),
        I0 = first(I),
        V0 = first(V)
      )
  )

  p_sim_test <- sim_test %>%
    pivot_longer(
      cols = c(T, L, I, V, C, CD4_total),
      names_to = "state",
      values_to = "value"
    ) %>%
    ggplot(aes(time, safe_log10(value), color = state)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(
      title = paste0(
        "Single-subject simulation test: id = ",
        first(sim_test$id)
      ),
      x = "Time",
      y = "log10(value)"
    )

  print(p_sim_test)

  save_plot_png_pdf(
    p_sim_test,
    "Fig0_single_subject_simulation_test",
    width = 8,
    height = 5
  )
}

# ============================================================
# 23. Optional debug M0 initial residuals
# ============================================================

flags_M0 <- c(
  za = 0,
  zk = 0,
  zxi = 0,
  zrho = 0
)

x0_M0 <- make_start_vector(
  model_flags = flags_M0,
  start_free_pars = default_start_free_pars
)

if (RUN_DEBUG_M0) {

  debug_joined_M0 <- debug_residuals(
    x = x0_M0,
    model_flags = flags_M0
  )

  if (!is.null(debug_joined_M0)) {
    write.csv(
      debug_joined_M0,
      file.path(out_dir, "debug_joined_M0_initial.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# 24. Optional test M0
# ============================================================

if (RUN_TEST_M0) {

  test_fit_M0 <- fit_one_model(
    model_row = candidate_models[1, ],
    start_free_pars = default_start_free_pars,
    start_id = 0
  )

  print(test_fit_M0$summary)
  print(test_fit_M0$params)

  write.csv(
    test_fit_M0$summary,
    file.path(out_dir, "test_fit_M0_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    test_fit_M0$params,
    file.path(out_dir, "test_fit_M0_params.csv"),
    row.names = FALSE
  )

  if (!is.null(test_fit_M0$fit)) {

    m0_residual_table <- make_residual_table(
      x = test_fit_M0$fit$par,
      model_flags = flags_M0
    )

    write.csv(
      m0_residual_table,
      file.path(out_dir, "test_fit_M0_residual_table.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# 25. Two-stage anchored random multi-start fitting
# ============================================================
# Stage 1:
#   Forward nested anchored fitting, same logic as the original script.
#   Parent models can provide anchors to child models.
#
# Stage 2:
#   Global reverse anchoring.
#   Best solutions found in Stage 1 from ALL models are converted to
#   starting values and passed back to EVERY model, including M0.
#
# Why this is necessary:
#   If a complex model has theta values close to zero but lower RSS than M0,
#   then that complex model is almost numerically equivalent to M0.
#   M0 should be given that parameter region as an anchor.
#   Otherwise AICc ranking can reflect optimizer luck rather than model support.

anchor_key <- function(start_vec, digits = 6) {

  positive_names <- c("a_L", "k_C", "xi_C", "rho_C")
  theta_names <- c("theta_a", "theta_k", "theta_xi", "theta_rho")

  positive_part <- round(log(pmax(as.numeric(start_vec[positive_names]), EPS)), digits)
  theta_part <- round(as.numeric(start_vec[theta_names]), digits)

  paste(
    c(positive_part, theta_part),
    collapse = "_"
  )
}

dedupe_anchor_starts <- function(anchor_starts) {

  if (length(anchor_starts) == 0) {
    return(list())
  }

  keys <- purrr::map_chr(anchor_starts, anchor_key)
  anchor_starts[!duplicated(keys)]
}

valid_start_free <- function(start_vec) {

  required_names <- c(
    "a_L",
    "k_C",
    "xi_C",
    "rho_C",
    "theta_a",
    "theta_k",
    "theta_xi",
    "theta_rho"
  )

  all(required_names %in% names(start_vec)) &&
    all(is.finite(start_vec)) &&
    all(start_vec[c("a_L", "k_C", "xi_C", "rho_C")] > 0)
}

collect_best_stage_anchors <- function(
    fit_list,
    max_anchors = Inf
) {

  anchor_tbl_list <- list()

  for (i in seq_along(fit_list)) {

    best_fit_i <- fit_list[[i]]$best_fit

    if (
      is.null(best_fit_i) ||
        is.null(best_fit_i$params) ||
        is.null(best_fit_i$summary) ||
        nrow(best_fit_i$params) == 0 ||
        nrow(best_fit_i$summary) == 0
    ) {
      next
    }

    param_row_i <- best_fit_i$params[1, ]
    summary_row_i <- best_fit_i$summary[1, ]

    if (
      is.na(param_row_i$a_L[1]) ||
        is.na(param_row_i$k_C[1]) ||
        is.na(param_row_i$xi_C[1]) ||
        is.na(param_row_i$rho_C[1]) ||
        is.na(summary_row_i$RSS[1]) ||
        !is.finite(summary_row_i$RSS[1])
    ) {
      next
    }

    start_vec_i <- params_row_to_start_free(param_row_i)

    if (!valid_start_free(start_vec_i)) {
      next
    }

    anchor_tbl_list[[length(anchor_tbl_list) + 1]] <- tibble(
      source_model = as.character(param_row_i$model[1]),
      source_RSS = as.numeric(summary_row_i$RSS[1]),
      source_AICc = as.numeric(summary_row_i$AICc[1]),
      a_L = as.numeric(start_vec_i["a_L"]),
      k_C = as.numeric(start_vec_i["k_C"]),
      xi_C = as.numeric(start_vec_i["xi_C"]),
      rho_C = as.numeric(start_vec_i["rho_C"]),
      theta_a = as.numeric(start_vec_i["theta_a"]),
      theta_k = as.numeric(start_vec_i["theta_k"]),
      theta_xi = as.numeric(start_vec_i["theta_xi"]),
      theta_rho = as.numeric(start_vec_i["theta_rho"]),
      anchor_key = anchor_key(start_vec_i)
    )
  }

  if (length(anchor_tbl_list) == 0) {
    return(list(
      anchor_table = tibble(),
      anchor_starts = list()
    ))
  }

  anchor_table <- bind_rows(anchor_tbl_list) %>%
    arrange(source_RSS, source_AICc) %>%
    distinct(anchor_key, .keep_all = TRUE)

  if (is.finite(max_anchors)) {
    anchor_table <- anchor_table %>%
      slice_head(n = max_anchors)
  }

  anchor_starts <- purrr::map(
    seq_len(nrow(anchor_table)),
    function(i) {
      c(
        a_L = anchor_table$a_L[i],
        k_C = anchor_table$k_C[i],
        xi_C = anchor_table$xi_C[i],
        rho_C = anchor_table$rho_C[i],
        theta_a = anchor_table$theta_a[i],
        theta_k = anchor_table$theta_k[i],
        theta_xi = anchor_table$theta_xi[i],
        theta_rho = anchor_table$theta_rho[i]
      )
    }
  )

  anchor_starts <- dedupe_anchor_starts(anchor_starts)

  list(
    anchor_table = anchor_table,
    anchor_starts = anchor_starts
  )
}

cat("\n============================================================\n")
cat("Stage 1 fitting: forward nested anchors + random starts\n")
cat("============================================================\n")

stage1_fit_list <- vector("list", nrow(candidate_models))

for (i in seq_len(nrow(candidate_models))) {

  model_row_i <- candidate_models[i, ]

  anchor_starts_i <- list()

  if (i > 1) {

    for (j in seq_len(i - 1)) {

      parent_row_j <- candidate_models[j, ]

      if (is_nested_parent(parent_row_j, model_row_i)) {

        parent_params_j <- stage1_fit_list[[j]]$best_fit$params

        if (
          nrow(parent_params_j) > 0 &&
            !is.na(parent_params_j$a_L[1]) &&
            !is.na(parent_params_j$k_C[1]) &&
            !is.na(parent_params_j$xi_C[1]) &&
            !is.na(parent_params_j$rho_C[1])
        ) {

          anchor_start_j <- params_row_to_start_free(parent_params_j[1, ])

          if (valid_start_free(anchor_start_j)) {
            anchor_starts_i[[length(anchor_starts_i) + 1]] <- anchor_start_j
          }
        }
      }
    }
  }

  anchor_starts_i <- dedupe_anchor_starts(anchor_starts_i)

  stage1_fit_list[[i]] <- fit_one_model_random_multistart(
    model_row = model_row_i,
    n_start = N_RANDOM_START_STAGE1,
    seed = RANDOM_SEED_BASE + i,
    anchor_starts = anchor_starts_i
  )
}

stage1_fit_list_best <- purrr::map(
  stage1_fit_list,
  "best_fit"
)

stage1_model_fit_table <- purrr::map_dfr(
  stage1_fit_list_best,
  "summary"
) %>%
  mutate(
    fitting_stage = "stage1_forward_nested"
  )

stage1_best_fit_params <- purrr::map_dfr(
  stage1_fit_list_best,
  "params"
) %>%
  mutate(
    fitting_stage = "stage1_forward_nested"
  )

stage1_all_random_fit_summaries <- purrr::map_dfr(
  seq_along(stage1_fit_list),
  function(i) {
    stage1_fit_list[[i]]$all_summaries
  }
) %>%
  mutate(
    fitting_stage = "stage1_forward_nested"
  )

stage2_anchor_object <- collect_best_stage_anchors(
  fit_list = stage1_fit_list,
  max_anchors = MAX_STAGE2_ANCHORS_PER_MODEL
)

stage2_global_anchor_table <- stage2_anchor_object$anchor_table
stage2_global_anchor_starts <- stage2_anchor_object$anchor_starts

write.csv(
  stage1_model_fit_table,
  file.path(out_dir, "stage1_model_fit_summary.csv"),
  row.names = FALSE
)

write.csv(
  stage1_best_fit_params,
  file.path(out_dir, "stage1_best_fit_params.csv"),
  row.names = FALSE
)

write.csv(
  stage1_all_random_fit_summaries,
  file.path(out_dir, "stage1_all_random_fit_summaries.csv"),
  row.names = FALSE
)

write.csv(
  stage2_global_anchor_table,
  file.path(out_dir, "stage2_global_anchor_table_from_stage1.csv"),
  row.names = FALSE
)

cat("\nStage 1 complete.\n")
cat("Number of global anchors collected for Stage 2:", length(stage2_global_anchor_starts), "\n")

if (RUN_TWO_STAGE_ANCHORED_FITTING) {

  cat("\n============================================================\n")
  cat("Stage 2 fitting: global reverse anchors + random starts\n")
  cat("============================================================\n")

  stage2_fit_list <- vector("list", nrow(candidate_models))

  for (i in seq_len(nrow(candidate_models))) {

    model_row_i <- candidate_models[i, ]

    # Every model receives the same Stage 1 global anchors.
    # make_start_vector() automatically ignores theta values not used by a target model.
    anchor_starts_i <- stage2_global_anchor_starts

    stage2_fit_list[[i]] <- fit_one_model_random_multistart(
      model_row = model_row_i,
      n_start = N_RANDOM_START_STAGE2,
      seed = RANDOM_SEED_BASE + 10000 + i,
      anchor_starts = anchor_starts_i
    )
  }

  stage2_fit_list_best <- purrr::map(
    stage2_fit_list,
    "best_fit"
  )

  stage2_model_fit_table <- purrr::map_dfr(
    stage2_fit_list_best,
    "summary"
  ) %>%
    mutate(
      fitting_stage = "stage2_global_reverse_anchor"
    )

  stage2_best_fit_params <- purrr::map_dfr(
    stage2_fit_list_best,
    "params"
  ) %>%
    mutate(
      fitting_stage = "stage2_global_reverse_anchor"
    )

  stage2_all_random_fit_summaries <- purrr::map_dfr(
    seq_along(stage2_fit_list),
    function(i) {
      stage2_fit_list[[i]]$all_summaries
    }
  ) %>%
    mutate(
      fitting_stage = "stage2_global_reverse_anchor"
    )

  write.csv(
    stage2_model_fit_table,
    file.path(out_dir, "stage2_model_fit_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    stage2_best_fit_params,
    file.path(out_dir, "stage2_best_fit_params.csv"),
    row.names = FALSE
  )

  write.csv(
    stage2_all_random_fit_summaries,
    file.path(out_dir, "stage2_all_random_fit_summaries.csv"),
    row.names = FALSE
  )

  # Final model-selection results are taken from Stage 2.
  random_fit_list <- stage2_fit_list

  fit_list_best <- stage2_fit_list_best

  model_fit_table_ms <- stage2_model_fit_table

  best_fit_params_ms <- stage2_best_fit_params

  all_random_fit_summaries <- bind_rows(
    stage1_all_random_fit_summaries,
    stage2_all_random_fit_summaries
  )

} else {

  # Fallback: use Stage 1 only.
  random_fit_list <- stage1_fit_list

  fit_list_best <- stage1_fit_list_best

  model_fit_table_ms <- stage1_model_fit_table

  best_fit_params_ms <- stage1_best_fit_params

  all_random_fit_summaries <- stage1_all_random_fit_summaries
}

write.csv(
  all_random_fit_summaries,
  file.path(out_dir, "all_random_fit_summaries_all_stages.csv"),
  row.names = FALSE
)

final_fitting_stage <- ifelse(
  RUN_TWO_STAGE_ANCHORED_FITTING,
  "stage2_global_reverse_anchor",
  "stage1_forward_nested"
)

all_random_fit_summaries_final_stage <- all_random_fit_summaries %>%
  filter(fitting_stage == final_fitting_stage)

# For backward compatibility with previous analysis scripts, this filename
# contains the final-stage summaries when two-stage fitting is enabled.
write.csv(
  all_random_fit_summaries_final_stage,
  file.path(out_dir, "all_random_fit_summaries.csv"),
  row.names = FALSE
)

cat("\nFinal fitting stage used for model selection:\n")
cat(final_fitting_stage, "\n")

# ============================================================
# 27. AICc ranking
# ============================================================

if (all(is.na(model_fit_table_ms$AICc))) {
  stop("All model fittings failed. Please check debug_residuals() output first.")
}

aicc_table_ms <- model_fit_table_ms %>%
  arrange(AICc) %>%
  mutate(
    delta_AICc = AICc - min(AICc, na.rm = TRUE),
    Akaike_weight = exp(-0.5 * delta_AICc) /
      sum(exp(-0.5 * delta_AICc), na.rm = TRUE),
    rank = row_number()
  )

# ============================================================
# 28. Nested RSS check
# ============================================================

nested_rss_check <- expand.grid(
  parent = candidate_models$model,
  child = candidate_models$model,
  stringsAsFactors = FALSE
) %>%
  left_join(
    candidate_models %>%
      rename(
        parent = model,
        parent_za = za,
        parent_zk = zk,
        parent_zxi = zxi,
        parent_zrho = zrho
      ),
    by = "parent"
  ) %>%
  left_join(
    candidate_models %>%
      rename(
        child = model,
        child_za = za,
        child_zk = zk,
        child_zxi = zxi,
        child_zrho = zrho
      ),
    by = "child"
  ) %>%
  filter(
    parent != child,
    parent_za <= child_za,
    parent_zk <= child_zk,
    parent_zxi <= child_zxi,
    parent_zrho <= child_zrho,
    parent_za < child_za |
      parent_zk < child_zk |
      parent_zxi < child_zxi |
      parent_zrho < child_zrho
  ) %>%
  left_join(
    model_fit_table_ms %>%
      select(model, parent_RSS = RSS),
    by = c("parent" = "model")
  ) %>%
  left_join(
    model_fit_table_ms %>%
      select(model, child_RSS = RSS),
    by = c("child" = "model")
  ) %>%
  mutate(
    child_minus_parent_RSS = child_RSS - parent_RSS,
    violation = child_minus_parent_RSS > 1e-6
  )

print(nested_rss_check)

if (any(nested_rss_check$violation, na.rm = TRUE)) {
  warning(
    "Some nested child models have larger RSS than parent models. ",
    "This indicates possible local minima or insufficient starts."
  )
}

# ============================================================
# 29. Most parsimonious model set
# ============================================================

most_parsimonious_models <- aicc_table_ms %>%
  filter(delta_AICc < 5) %>%
  arrange(AICc)

# ============================================================
# 30. Parameter boundary check
# ============================================================

boundary_check <- best_fit_params_ms %>%
  mutate(
    a_L_near_lower = !is.na(a_L) & a_L <= 1.01e-8,
    a_L_near_upper = !is.na(a_L) & a_L >= 0.99,
    k_C_near_lower = !is.na(k_C) & k_C <= 1.01e-6,
    k_C_near_upper = !is.na(k_C) & k_C >= 9.9,
    xi_C_near_lower = !is.na(xi_C) & xi_C <= 1.01e-5,
    xi_C_near_upper = !is.na(xi_C) & xi_C >= 99,
    rho_C_near_lower = !is.na(rho_C) & rho_C <= 1.01e-6,
    rho_C_near_upper = !is.na(rho_C) & rho_C >= 4.95,
    theta_a_near_boundary = !is.na(theta_a) & abs(theta_a) >= 4.95,
    theta_k_near_boundary = !is.na(theta_k) & abs(theta_k) >= 4.95,
    theta_xi_near_boundary = !is.na(theta_xi) & abs(theta_xi) >= 4.95,
    theta_rho_near_boundary = !is.na(theta_rho) & abs(theta_rho) >= 4.95
  )

print(boundary_check)

if (
  any(
    boundary_check$a_L_near_lower |
      boundary_check$a_L_near_upper |
      boundary_check$k_C_near_lower |
      boundary_check$k_C_near_upper |
      boundary_check$xi_C_near_lower |
      boundary_check$xi_C_near_upper |
      boundary_check$rho_C_near_lower |
      boundary_check$rho_C_near_upper |
      boundary_check$theta_a_near_boundary |
      boundary_check$theta_k_near_boundary |
      boundary_check$theta_xi_near_boundary |
      boundary_check$theta_rho_near_boundary,
    na.rm = TRUE
  )
) {
  warning(
    "Some fitted parameters are close to bounds. ",
    "Interpret these models with caution."
  )
}

# ============================================================
# 31. Print main results
# ============================================================

print(best_fit_params_ms)

aicc_table_ms %>%
  select(
    rank,
    model,
    RSS,
    raw_RSS,
    logLik,
    AICc,
    delta_AICc,
    Akaike_weight,
    q,
    K,
    best_start_id,
    n_random_start,
    n_anchor_start,
    convergence_info
  ) %>%
  print(n = Inf)

most_parsimonious_models %>%
  select(
    rank,
    model,
    RSS,
    raw_RSS,
    logLik,
    AICc,
    delta_AICc,
    Akaike_weight,
    q,
    K
  ) %>%
  print(n = Inf)

# ============================================================
# 32. Save model selection results
# ============================================================

write.csv(
  best_fit_params_ms,
  file.path(out_dir, "best_fit_params_multistart.csv"),
  row.names = FALSE
)

write.csv(
  model_fit_table_ms,
  file.path(out_dir, "model_fit_summary_multistart.csv"),
  row.names = FALSE
)

write.csv(
  aicc_table_ms,
  file.path(out_dir, "aicc_table_multistart.csv"),
  row.names = FALSE
)

write.csv(
  most_parsimonious_models,
  file.path(out_dir, "most_parsimonious_models_deltaAICc_lt5.csv"),
  row.names = FALSE
)

write.csv(
  all_random_fit_summaries_final_stage,
  file.path(out_dir, "all_random_fit_summaries.csv"),
  row.names = FALSE
)

write.csv(
  nested_rss_check,
  file.path(out_dir, "nested_rss_check.csv"),
  row.names = FALSE
)

write.csv(
  boundary_check,
  file.path(out_dir, "parameter_boundary_check.csv"),
  row.names = FALSE
)

saveRDS(
  random_fit_list,
  file.path(out_dir, "random_fit_list.rds")
)

# ============================================================
# 33. Best model residuals, predictions, simulations
# ============================================================

best_model_name <- aicc_table_ms %>%
  arrange(AICc) %>%
  slice(1) %>%
  pull(model)

best_model_index <- match(
  best_model_name,
  candidate_models$model
)

best_fit_overall <- fit_list_best[[best_model_index]]

if (!is.null(best_fit_overall$fit)) {

  best_model_row <- candidate_models[best_model_index, ]

  best_model_flags <- c(
    za = as.numeric(best_model_row$za),
    zk = as.numeric(best_model_row$zk),
    zxi = as.numeric(best_model_row$zxi),
    zrho = as.numeric(best_model_row$zrho)
  )

  best_x <- best_fit_overall$fit$par

  best_residual_table <- make_residual_table(
    x = best_x,
    model_flags = best_model_flags
  )

  write.csv(
    best_residual_table,
    file.path(out_dir, "best_model_residual_table.csv"),
    row.names = FALSE
  )

  best_pars <- unpack_pars(
    x = best_x,
    model_flags = best_model_flags
  )

  best_sim <- simulate_all(
    pars = best_pars,
    model_flags = best_model_flags
  )

  best_pred <- make_predictions(best_sim)

  write.csv(
    best_sim,
    file.path(out_dir, "best_model_simulation_full.csv"),
    row.names = FALSE
  )

  write.csv(
    best_pred,
    file.path(out_dir, "best_model_predictions_long.csv"),
    row.names = FALSE
  )

  cd8_effect_summary <- best_sim %>%
    group_by(id, group) %>%
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
      a_L_G = median(a_L_G, na.rm = TRUE),
      p_G = median(p_G, na.rm = TRUE),
      k_C_G = median(k_C_G, na.rm = TRUE),
      xi_C_G = median(xi_C_G, na.rm = TRUE),
      rho_C_G = median(rho_C_G, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(group, id)

  write.csv(
    cd8_effect_summary,
    file.path(out_dir, "best_model_CD8_effect_summary_by_id.csv"),
    row.names = FALSE
  )

  group_bias_summary <- best_residual_table %>%
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
    group_bias_summary,
    file.path(out_dir, "best_model_group_bias_summary.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 34. Save run settings
# ============================================================

run_settings <- tibble(
  item = c(
    "N_RANDOM_START_STAGE1",
    "N_RANDOM_START_STAGE2",
    "RUN_TWO_STAGE_ANCHORED_FITTING",
    "MAX_STAGE2_ANCHORS_PER_MODEL",
    "final_fitting_stage",
    "RANDOM_SEED_BASE",
    "USE_VARIABLE_SCALING",
    "FIT_VARIABLES",
    "data_file",
    "init_file",
    "fixed_p",
    "fixed_K_C",
    "fixed_K_kill",
    "CD8_scale",
    "fixed_d_T",
    "fixed_beta",
    "fixed_f_L",
    "fixed_d_L",
    "fixed_delta0",
    "fixed_c",
    "fixed_d_C",
    "theta_sign_convention",
    "model_structure",
    "model_selection_parameters",
    "CD4_handling",
    "CD8_handling",
    "DNA_handling",
    "pVL_LOD_handling",
    "out_dir"
  ),
  value = c(
    as.character(N_RANDOM_START_STAGE1),
    as.character(N_RANDOM_START_STAGE2),
    as.character(RUN_TWO_STAGE_ANCHORED_FITTING),
    as.character(MAX_STAGE2_ANCHORS_PER_MODEL),
    final_fitting_stage,
    as.character(RANDOM_SEED_BASE),
    as.character(USE_VARIABLE_SCALING),
    paste(FIT_VARIABLES, collapse = ","),
    DATA_FILE,
    INIT_FILE,
    as.character(fixed_pars["p"]),
    as.character(fixed_pars["K_C"]),
    as.character(fixed_pars["K_kill"]),
    as.character(fixed_pars["CD8_scale"]),
    as.character(fixed_pars["d_T"]),
    as.character(fixed_pars["beta"]),
    as.character(fixed_pars["f_L"]),
    as.character(fixed_pars["d_L"]),
    as.character(fixed_pars["delta0"]),
    as.character(fixed_pars["c"]),
    as.character(fixed_pars["d_C"]),
    "parameter_early = parameter_late * exp(theta); theta > 0 means early > late",
    "T,L,I,V,C; no exhaustion compartment; saturated CD8 killing; fixed p; new C_base_group CD8 equation",
    "a_L, k_C, xi_C, rho_C",
    "CD4 raw cells/mL converted to log10; model prediction log10(T + L + I)",
    "CD8 raw cells/mL converted to log10; model prediction log10(C); dC uses C_base_group, not C0",
    "DNA already log10; model prediction log10(1e6 * (L + I)/(T + L + I))",
    "pVL value and LOD are raw copies/mL; both converted to log10; censoring judged on raw scale",
    out_dir
  )
)

write.csv(
  run_settings,
  file.path(out_dir, "run_settings.csv"),
  row.names = FALSE
)

# ============================================================
# 35. Figures
# ============================================================

p_aicc <- aicc_table_ms %>%
  mutate(
    model = factor(model, levels = rev(model))
  ) %>%
  ggplot(aes(x = model, y = delta_AICc)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Model selection by AICc: fixed p = 8000",
    x = "Model",
    y = expression(Delta * "AICc")
  ) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  annotate(
    "text",
    x = 1,
    y = 5,
    label = expression(Delta * "AICc = 5"),
    hjust = -0.05,
    vjust = -0.5,
    size = 3.5
  )

print(p_aicc)

save_plot_png_pdf(
  p_aicc,
  "Fig1_AICc_model_ranking",
  width = 8,
  height = 6
)

if (exists("best_residual_table") && !is.null(best_residual_table)) {

  p_obs_pred <- best_residual_table %>%
    ggplot(aes(x = value, y = pred)) +
    geom_point(aes(shape = is_censored), alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ variable, scales = "free") +
    theme_bw() +
    labs(
      title = paste0("Observed vs predicted values: ", best_model_name),
      x = "Observed value used for fitting",
      y = "Predicted value",
      shape = "Censored"
    )

  print(p_obs_pred)

  save_plot_png_pdf(
    p_obs_pred,
    "Fig2_observed_vs_predicted_best_model",
    width = 9,
    height = 5
  )

  p_resid <- best_residual_table %>%
    ggplot(aes(x = variable, y = resid)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +
    theme_bw() +
    labs(
      title = paste0("Weighted residuals by variable: ", best_model_name),
      x = "Variable",
      y = "Weighted residual"
    )

  print(p_resid)

  save_plot_png_pdf(
    p_resid,
    "Fig3_residual_by_variable_best_model",
    width = 7,
    height = 5
  )

  p_timecourse <- best_residual_table %>%
    ggplot(aes(x = time)) +
    geom_line(
      aes(y = pred, group = id),
      alpha = 0.35,
      linewidth = 0.5
    ) +
    geom_point(
      aes(y = value, shape = is_censored),
      alpha = 0.7,
      size = 1.8
    ) +
    facet_grid(variable ~ group, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Observed and predicted time courses: ", best_model_name),
      x = "Time",
      y = "Value used for fitting",
      shape = "Censored"
    )

  print(p_timecourse)

  save_plot_png_pdf(
    p_timecourse,
    "Fig4_timecourse_observed_predicted_best_model",
    width = 10,
    height = 8
  )

  rss_by_variable <- best_residual_table %>%
    group_by(variable) %>%
    summarise(
      n = n(),
      n_censored = sum(is_censored, na.rm = TRUE),
      raw_RSS = sum(raw_resid^2, na.rm = TRUE),
      weighted_RSS = sum(resid^2, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(
    rss_by_variable,
    file.path(out_dir, "best_model_RSS_by_variable.csv"),
    row.names = FALSE
  )

  p_rss <- rss_by_variable %>%
    ggplot(aes(x = variable, y = weighted_RSS)) +
    geom_col(width = 0.65) +
    theme_bw() +
    labs(
      title = paste0("Weighted RSS contribution by variable: ", best_model_name),
      x = "Variable",
      y = "Weighted RSS"
    )

  print(p_rss)

  save_plot_png_pdf(
    p_rss,
    "Fig5_weighted_RSS_by_variable",
    width = 7,
    height = 5
  )
}

if (exists("best_fit_overall") && !is.null(best_fit_overall$params)) {

  best_param_plot_data <- best_fit_overall$params %>%
    select(
      model,
      a_L_late,
      a_L_early,
      k_C_late,
      k_C_early,
      xi_C_late,
      xi_C_early,
      rho_C_late,
      rho_C_early
    ) %>%
    pivot_longer(
      cols = -model,
      names_to = "parameter_group",
      values_to = "value"
    ) %>%
    mutate(
      parameter = case_when(
        grepl("^a_L", parameter_group) ~ "a_L",
        grepl("^k_C", parameter_group) ~ "k_C",
        grepl("^xi_C", parameter_group) ~ "xi_C",
        grepl("^rho_C", parameter_group) ~ "rho_C",
        TRUE ~ parameter_group
      ),
      group = case_when(
        grepl("_early$", parameter_group) ~ "early",
        grepl("_late$", parameter_group) ~ "late",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))

  p_param <- best_param_plot_data %>%
    ggplot(aes(x = group, y = value)) +
    geom_col(width = 0.65) +
    facet_wrap(~ parameter, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Early vs late parameters: ", best_model_name),
      x = "Group",
      y = "Estimated parameter value"
    )

  print(p_param)

  save_plot_png_pdf(
    p_param,
    "Fig6_best_model_parameter_comparison",
    width = 10,
    height = 5
  )
}

if (exists("best_sim") && !is.null(best_sim)) {

  cd8_timecourse_plot_data <- best_sim %>%
    select(
      id,
      group,
      time,
      C,
      C_base_group,
      C_scaled,
      CD8_activation_fraction,
      CD8_killing_saturation,
      CD8_killing_rate,
      p_eff_over_p_G
    ) %>%
    pivot_longer(
      cols = c(
        C,
        C_base_group,
        C_scaled,
        CD8_activation_fraction,
        CD8_killing_saturation,
        CD8_killing_rate,
        p_eff_over_p_G
      ),
      names_to = "quantity",
      values_to = "value"
    )

  p_cd8_time <- cd8_timecourse_plot_data %>%
    ggplot(aes(x = time, y = value, group = id)) +
    geom_line(alpha = 0.4, linewidth = 0.5) +
    facet_grid(quantity ~ group, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("CD8 state and effective functions: ", best_model_name),
      x = "Time",
      y = "Value"
    )

  print(p_cd8_time)

  save_plot_png_pdf(
    p_cd8_time,
    "Fig7_CD8_states_and_effects_timecourse",
    width = 11,
    height = 9
  )
}

if (exists("best_residual_table") && !is.null(best_residual_table)) {

  pVL_residual_by_group_and_censoring <- best_residual_table %>%
    filter(variable == "pVL") %>%
    group_by(group, is_censored) %>%
    summarise(
      n = n(),
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
    pVL_residual_by_group_and_censoring,
    file.path(out_dir, "pVL_residual_by_group_and_censoring.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 36. Final message
# ============================================================

cat("\nAll results saved to:\n")
cat(out_dir, "\n")

cat("\nFigures saved to:\n")
cat(fig_dir, "\n")

cat("\nMain files to check:\n")
cat("1. candidate_models.csv\n")
cat("2. init_table_with_C0_and_Cbase_used.csv\n")
cat("3. data_scale_check.csv\n")
cat("4. aicc_table_multistart.csv\n")
cat("5. best_fit_params_multistart.csv\n")
cat("5a. stage1_model_fit_summary.csv\n")
cat("5b. stage2_model_fit_summary.csv\n")
cat("5c. stage2_global_anchor_table_from_stage1.csv\n")
cat("5d. all_random_fit_summaries_all_stages.csv\n")
cat("6. nested_rss_check.csv\n")
cat("7. parameter_boundary_check.csv\n")
cat("8. best_model_residual_table.csv\n")
cat("9. best_model_RSS_by_variable.csv\n")
cat("10. best_model_CD8_effect_summary_by_id.csv\n")
cat("11. best_model_group_bias_summary.csv\n")
cat("12. pVL_residual_by_group_and_censoring.csv\n")
cat("13. run_settings.csv\n")

cat("\nMain figures to show:\n")
cat("1. Fig1_AICc_model_ranking.png\n")
cat("2. Fig2_observed_vs_predicted_best_model.png\n")
cat("3. Fig3_residual_by_variable_best_model.png\n")
cat("4. Fig4_timecourse_observed_predicted_best_model.png\n")
cat("5. Fig5_weighted_RSS_by_variable.png\n")
cat("6. Fig6_best_model_parameter_comparison.png\n")
cat("7. Fig7_CD8_states_and_effects_timecourse.png\n")

cat("\nBest model:\n")
cat(best_model_name, "\n")
