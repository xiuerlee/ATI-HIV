library(deSolve)
library(tidyverse)
library(minpack.lm)

source("model_with_CD8count.R")


# ============================================================
# 0. 全局设置
# ============================================================

EPS <- 1e-12

N_RANDOM_START <- 5
RANDOM_SEED_BASE <- 1000

USE_VARIABLE_SCALING <- TRUE

# 当前 noDNA 版本只使用 pVL 和 CD4 进入残差、RSS、AICc 和模型选择。
# DNA 不作为拟合目标。
FIT_VARIABLES <- c("pVL", "CD4")

RUN_SINGLE_SIM_TEST <- TRUE
RUN_DEBUG_M0 <- TRUE
RUN_TEST_M0 <- TRUE

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  "../Rdata/model_selection_results_CD8count_noDNA",
  run_tag
)

fig_dir <- file.path(out_dir, "figures")

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
# 1. 读入数据
# ============================================================
# 注意：
# 当前代码默认 raw_data_with_cd8_counts.csv 中：
# CD4 已经是 cells/mL；
# CD8 已经是 cells/mL；
# pVL 是 raw copies/mL；
# pVL 的 lod 也是 raw copies/mL；
# DNA 即使存在于数据表中，本 noDNA 版本也不用于拟合。

data_clean <- read.csv("../Rdata/raw_data_with_cd8_counts.csv")
init_table <- read.csv("../Rdata/init_table_I0_replaced.csv")


# ============================================================
# 2. 固定参数
# ============================================================
# delta0 固定。
# CD8 杀伤写成 k_CD8 * CD8_count_G * I。

fixed_pars <- c(
  d_T = 0.01,
  beta = 1.58e-8,
  f_L = 0.005,
  d_L = 5e-4,
  delta0 = 1.4,
  c = 23
)


# ============================================================
# 3. 初始拟合参数
# ============================================================
# 注意：
# 因为 CD8_count_G 现在是 cells/mL，
# k_CD8 必须比 cells/uL 尺度下小 1000 倍左右。
#
# 若 CD8_count_G ≈ 1e6 cells/mL，
# k_CD8 = 1e-8 对应有效 killing rate ≈ 0.01 day^-1。

default_start_free_pars <- c(
  a_L = 1e-3,
  p = 6e3,
  k_CD8 = 1e-8,
  theta_a = 0,
  theta_p = 0,
  theta_k = 0
)


# ============================================================
# 4. 测试用参数
# ============================================================

test_pars <- c(
  fixed_pars,
  a_L = 1e-3,
  p = 6e3,
  k_CD8 = 1e-8,
  theta_a = 0,
  theta_p = 0,
  theta_k = 0
)

test_flags <- c(
  za = 0,
  zp = 0,
  zk = 0
)


# ============================================================
# 5. 工具函数：组别识别
# ============================================================

is_early_group <- function(group_value) {
  tolower(as.character(group_value)) %in% c(
    "early",
    "early_art",
    "early art",
    "early-treated",
    "early treated",
    "et"
  )
}


# ============================================================
# 6. 工具函数：保存图
# ============================================================

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


# ============================================================
# 7. 数据基础清理
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
    V0 = as.numeric(V0)
  ) %>%
  filter(
    !is.na(T0),
    !is.na(L0),
    !is.na(V0)
  )


# ============================================================
# 8. CD4 / CD8 单位检查
# ============================================================
# 这里不修改数据，只做检查。
# 如果 CD4 / CD8 中位数仍然是几百到几千，
# 说明它们可能还是 cells/uL，不是 cells/mL。

unit_check <- data_raw %>%
  filter(variable %in% c(
    "CD4",
    "CD8",
    "CD8_count",
    "CD8_counts",
    "cd8",
    "cd8_count",
    "cd8_counts"
  )) %>%
  group_by(variable) %>%
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
  file.path(out_dir, "CD4_CD8_unit_check.csv"),
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
    unit_check$variable %in% c("CD8", "CD8_count", "CD8_counts", "cd8", "cd8_count", "cd8_counts") &
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
# 9. 提取 CD8 count
# ============================================================
# CD8 是外源协变量，不作为拟合目标。
# 优先使用每个个体 time == 0 的 CD8。
# 如果没有 time == 0，则使用该个体所有 CD8 的中位数。
# 如果该个体没有 CD8，则退回到同组中位数。
# 如果同组也没有，则退回到全体中位数。
#
# 注意：
# 这里默认 CD8 的 value 已经是 cells/mL，不再乘以 1000。

cd8_variable_candidates <- c(
  "CD8",
  "CD8_count",
  "CD8_counts",
  "cd8",
  "cd8_count",
  "cd8_counts"
)

cd8_data <- data_raw %>%
  filter(
    variable %in% cd8_variable_candidates,
    !is.na(value),
    is.finite(value)
  )

if (nrow(cd8_data) == 0) {
  stop(
    "No CD8 count data found. Please check whether the CD8 variable is named as one of: ",
    paste(cd8_variable_candidates, collapse = ", ")
  )
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
    CD8_count_id = ifelse(
      is.finite(CD8_time0),
      CD8_time0,
      CD8_median_id
    )
  ) %>%
  select(id, group, CD8_count_id)

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
    CD8_count_G = case_when(
      is.finite(CD8_count_id) ~ CD8_count_id,
      is.finite(CD8_count_group) ~ CD8_count_group,
      TRUE ~ overall_cd8_median
    )
  ) %>%
  filter(
    is.finite(CD8_count_G),
    CD8_count_G > 0
  )

valid_ids <- init_fit$id


# ============================================================
# 10. 强制指定 pVL / CD4 / CD8 的尺度
# ============================================================
# pVL:
# value 和 lod 都是 raw copies/mL，都需要 log10 转换。
#
# CD4:
# 当前表格里的 CD4 视为已经是 cells/mL，作为 raw count，再 log10。
#
# CD8:
# 当前表格里的 CD8 视为已经是 cells/mL，作为外源协变量，不进入残差。
#
# DNA:
# 当前 noDNA 版本不拟合 DNA。即使 data_raw 中存在 DNA，也不会进入 data_fit、
# residual_scale_table、RSS、AICc 或最佳模型拟合图。

scale_settings <- tibble(
  variable = c("pVL", "CD4", "CD8", "DNA"),
  value_scale_used = c(
    "raw copies/mL -> log10",
    "raw cells/mL -> log10",
    "raw cells/mL used as external covariate",
    "excluded from fitting in this noDNA version"
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


# ============================================================
# 11. 拟合数据预处理
# ============================================================
# 关键修正：
# 1. pVL value_raw 和 lod_raw 都是 raw copies/mL，二者都要 log10。
# 2. pVL censored 判断使用 raw scale: value_raw <= lod_raw。
# 3. CD4 value_raw 已经是 cells/mL，因此直接 log10(value_raw)。
# 4. DNA 不进入拟合数据。

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
      variable == "pVL" ~ log10(pmax(value_raw, EPS)),
      variable == "CD4" ~ log10(pmax(value_raw, EPS)),
      TRUE ~ value_raw
    ),
    
    lod = case_when(
      variable == "pVL" &
        !is.na(lod_raw) &
        is.finite(lod_raw) ~ log10(pmax(lod_raw, EPS)),
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


# ============================================================
# 12. 检查数据尺度
# ============================================================

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
# 13. 残差尺度标准化
# ============================================================
# 目的：
# 防止某一个变量因为数值尺度更大而主导 RSS。
# AICc 使用的是标准化后的 RSS。

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
# 14. 8 个候选模型
# ============================================================
# za: 是否允许 a_L 有 early / late 差异。
# zp: 是否允许 p 有 early / late 差异。
# zk: 是否允许 k_CD8 有 early / late 差异。

candidate_models <- tibble::tribble(
  ~model,      ~za, ~zp, ~zk,
  "M0",          0,   0,   0,
  "Ma",          1,   0,   0,
  "Mp",          0,   1,   0,
  "Mk",          0,   0,   1,
  "Ma_p",        1,   1,   0,
  "Ma_k",        1,   0,   1,
  "Mp_k",        0,   1,   1,
  "Ma_p_k",      1,   1,   1
)


# ============================================================
# 15. 单个个体模拟函数
# ============================================================

simulate_one <- function(id_row, pars, times, model_flags) {
  
  G <- ifelse(is_early_group(id_row$group), 1, 0)
  
  a_L_G <- as.numeric(pars["a_L"]) *
    exp(-as.numeric(model_flags["za"]) * as.numeric(pars["theta_a"]) * G)
  
  p_G <- as.numeric(pars["p"]) *
    exp(-as.numeric(model_flags["zp"]) * as.numeric(pars["theta_p"]) * G)
  
  k_CD8_G <- as.numeric(pars["k_CD8"]) *
    exp(as.numeric(model_flags["zk"]) * as.numeric(pars["theta_k"]) * G)
  
  V0 <- max(as.numeric(id_row$V0), EPS)
  L0 <- max(as.numeric(id_row$L0), EPS)
  T0 <- max(as.numeric(id_row$T0), EPS)
  
  CD8_count_G <- max(as.numeric(id_row$CD8_count_G), EPS)
  
  I0 <- as.numeric(pars["c"]) * V0 / p_G
  I0 <- max(I0, EPS)
  
  state <- c(
    T = T0,
    L = L0,
    I = I0,
    V = V0
  )
  
  pars_use <- pars
  pars_use["k_CD8"] <- k_CD8_G
  
  pars_use <- c(
    pars_use,
    T0 = T0,
    a_L_G = a_L_G,
    p_G = p_G,
    CD8_count_G = CD8_count_G
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
    message("ODE error for id ", id_row$id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) {
    return(NULL)
  }
  
  as.data.frame(out) %>%
    mutate(
      id = id_row$id,
      group = id_row$group,
      CD8_count_G = CD8_count_G,
      a_L_G = a_L_G,
      p_G = p_G,
      k_CD8_G = k_CD8_G,
      CD8_killing_rate_G = k_CD8_G * CD8_count_G,
      CD4_total = T + L + I,
      pred_pVL_log10 = log10(pmax(V, EPS)),
      pred_DNA_log10 = log10(pmax(1e6 * (L + I) / CD4_total, EPS)),
      pred_CD4_log10 = log10(pmax(CD4_total, EPS))
    )
}


# ============================================================
# 16. 根据候选模型生成自由参数初值
# ============================================================

make_start_vector <- function(model_flags, start_free_pars = default_start_free_pars) {
  
  x0 <- c(
    log_a_L = log(unname(start_free_pars["a_L"])),
    log_p = log(unname(start_free_pars["p"])),
    log_k_CD8 = log(unname(start_free_pars["k_CD8"]))
  )
  
  if (model_flags["za"] == 1) {
    x0 <- c(x0, theta_a = unname(start_free_pars["theta_a"]))
  }
  
  if (model_flags["zp"] == 1) {
    x0 <- c(x0, theta_p = unname(start_free_pars["theta_p"]))
  }
  
  if (model_flags["zk"] == 1) {
    x0 <- c(x0, theta_k = unname(start_free_pars["theta_k"]))
  }
  
  x0
}


# ============================================================
# 17. 参数下界
# ============================================================

make_lower_vector <- function(model_flags) {
  
  lower <- c(
    log_a_L = log(1e-8),
    log_p = log(1),
    log_k_CD8 = log(1e-12)
  )
  
  if (model_flags["za"] == 1) {
    lower <- c(lower, theta_a = -5)
  }
  
  if (model_flags["zp"] == 1) {
    lower <- c(lower, theta_p = -5)
  }
  
  if (model_flags["zk"] == 1) {
    lower <- c(lower, theta_k = -5)
  }
  
  lower
}


# ============================================================
# 18. 参数上界
# ============================================================

make_upper_vector <- function(model_flags) {
  
  upper <- c(
    log_a_L = log(1),
    log_p = log(1e8),
    log_k_CD8 = log(1e-5)
  )
  
  if (model_flags["za"] == 1) {
    upper <- c(upper, theta_a = 5)
  }
  
  if (model_flags["zp"] == 1) {
    upper <- c(upper, theta_p = 5)
  }
  
  if (model_flags["zk"] == 1) {
    upper <- c(upper, theta_k = 5)
  }
  
  upper
}


# ============================================================
# 19. 把 nls.lm 的自由参数转回模型参数
# ============================================================

unpack_pars <- function(x, model_flags) {
  
  pars <- c(
    fixed_pars,
    a_L = exp(unname(x["log_a_L"])),
    p = exp(unname(x["log_p"])),
    k_CD8 = exp(unname(x["log_k_CD8"])),
    theta_a = 0,
    theta_p = 0,
    theta_k = 0
  )
  
  if (model_flags["za"] == 1) {
    pars["theta_a"] <- unname(x["theta_a"])
  }
  
  if (model_flags["zp"] == 1) {
    pars["theta_p"] <- unname(x["theta_p"])
  }
  
  if (model_flags["zk"] == 1) {
    pars["theta_k"] <- unname(x["theta_k"])
  }
  
  pars
}


# ============================================================
# 20. 模拟全部个体
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
# 21. 预测函数
# ============================================================

make_predictions <- function(sim) {
  
  sim %>%
    select(
      id,
      group,
      time,
      pred_pVL_log10,
      pred_CD4_log10
    ) %>%
    pivot_longer(
      cols = starts_with("pred_"),
      names_to = "variable",
      values_to = "pred"
    ) %>%
    mutate(
      variable = recode(
        variable,
        pred_pVL_log10 = "pVL",
        pred_CD4_log10 = "CD4"
      )
    )
}


# ============================================================
# 22. 构建残差表
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


# ============================================================
# 23. 残差函数
# ============================================================

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
# 24. 诊断残差函数
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
        min_CD8 = min(CD8_count_G, na.rm = TRUE),
        max_CD8 = max(CD8_count_G, na.rm = TRUE),
        min_CD8_killing_rate = min(CD8_killing_rate_G, na.rm = TRUE),
        max_CD8_killing_rate = max(CD8_killing_rate_G, na.rm = TRUE),
        any_nonfinite = any(
          !is.finite(T) |
            !is.finite(L) |
            !is.finite(I) |
            !is.finite(V)
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
# 25. 拟合单个候选模型
# ============================================================

fit_one_model <- function(
    model_row,
    start_free_pars = default_start_free_pars,
    start_id = NA_integer_
) {
  
  model_flags <- c(
    za = as.numeric(model_row$za),
    zp = as.numeric(model_row$zp),
    zk = as.numeric(model_row$zk)
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
  
  if (is.null(fit)) {
    
    summary_row <- tibble(
      model = model_row$model,
      za = model_row$za,
      zp = model_row$zp,
      zk = model_row$zk,
      n = nrow(data_fit),
      q = length(x0),
      K = length(x0) + 1,
      RSS = NA_real_,
      raw_RSS = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      convergence_info = NA_integer_,
      start_id = start_id
    )
    
    param_row <- tibble(
      model = model_row$model,
      a_L = NA_real_,
      p = NA_real_,
      k_CD8 = NA_real_,
      theta_a = NA_real_,
      theta_p = NA_real_,
      theta_k = NA_real_,
      a_L_late = NA_real_,
      a_L_early = NA_real_,
      p_late = NA_real_,
      p_early = NA_real_,
      k_CD8_late = NA_real_,
      k_CD8_early = NA_real_,
      start_id = start_id
    )
    
    return(list(
      fit = NULL,
      summary = summary_row,
      params = param_row
    ))
  }
  
  x_hat <- fit$par
  pars_hat <- unpack_pars(x_hat, model_flags)
  
  residual_table_hat <- make_residual_table(
    x = x_hat,
    model_flags = model_flags
  )
  
  if (is.null(residual_table_hat)) {
    
    summary_row <- tibble(
      model = model_row$model,
      za = model_row$za,
      zp = model_row$zp,
      zk = model_row$zk,
      n = nrow(data_fit),
      q = length(x0),
      K = length(x0) + 1,
      RSS = NA_real_,
      raw_RSS = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      convergence_info = fit$info,
      start_id = start_id
    )
    
    param_row <- tibble(
      model = model_row$model,
      a_L = NA_real_,
      p = NA_real_,
      k_CD8 = NA_real_,
      theta_a = NA_real_,
      theta_p = NA_real_,
      theta_k = NA_real_,
      a_L_late = NA_real_,
      a_L_early = NA_real_,
      p_late = NA_real_,
      p_early = NA_real_,
      k_CD8_late = NA_real_,
      k_CD8_early = NA_real_,
      start_id = start_id
    )
    
    return(list(
      fit = fit,
      summary = summary_row,
      params = param_row
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
    exp(-unname(model_flags["za"]) * unname(pars_hat["theta_a"]))
  
  p_late <- unname(pars_hat["p"])
  p_early <- unname(pars_hat["p"]) *
    exp(-unname(model_flags["zp"]) * unname(pars_hat["theta_p"]))
  
  k_CD8_late <- unname(pars_hat["k_CD8"])
  k_CD8_early <- unname(pars_hat["k_CD8"]) *
    exp(unname(model_flags["zk"]) * unname(pars_hat["theta_k"]))
  
  summary_row <- tibble(
    model = model_row$model,
    za = model_row$za,
    zp = model_row$zp,
    zk = model_row$zk,
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
    p = unname(pars_hat["p"]),
    k_CD8 = unname(pars_hat["k_CD8"]),
    theta_a = unname(pars_hat["theta_a"]),
    theta_p = unname(pars_hat["theta_p"]),
    theta_k = unname(pars_hat["theta_k"]),
    a_L_late = a_L_late,
    a_L_early = a_L_early,
    p_late = p_late,
    p_early = p_early,
    k_CD8_late = k_CD8_late,
    k_CD8_early = k_CD8_early,
    start_id = start_id
  )
  
  list(
    fit = fit,
    summary = summary_row,
    params = param_row
  )
}


# ============================================================
# 26. 随机初值生成函数
# ============================================================
# 注意：
# k_CD8 是 cells/mL 尺度下的 killing coefficient。
# 随机范围设为 1e-11 到 1e-6。

make_random_start <- function() {
  
  c(
    a_L = exp(runif(1, log(1e-5), log(1e-1))),
    p = exp(runif(1, log(1e2), log(1e5))),
    k_CD8 = exp(runif(1, log(1e-11), log(1e-6))),
    theta_a = runif(1, -2, 2),
    theta_p = runif(1, -2, 2),
    theta_k = runif(1, -2, 2)
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


# ============================================================
# 27. anchored multi-start fitting
# ============================================================

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
# 28. 父模型 best-fit 转为子模型 anchor 初值
# ============================================================

params_row_to_start_free <- function(param_row) {
  
  c(
    a_L = as.numeric(param_row$a_L),
    p = as.numeric(param_row$p),
    k_CD8 = as.numeric(param_row$k_CD8),
    theta_a = ifelse(is.na(param_row$theta_a), 0, as.numeric(param_row$theta_a)),
    theta_p = ifelse(is.na(param_row$theta_p), 0, as.numeric(param_row$theta_p)),
    theta_k = ifelse(is.na(param_row$theta_k), 0, as.numeric(param_row$theta_k))
  )
}


is_nested_parent <- function(parent_row, child_row) {
  
  parent_flags <- c(
    za = as.numeric(parent_row$za),
    zp = as.numeric(parent_row$zp),
    zk = as.numeric(parent_row$zk)
  )
  
  child_flags <- c(
    za = as.numeric(child_row$za),
    zp = as.numeric(child_row$zp),
    zk = as.numeric(child_row$zk)
  )
  
  all(parent_flags <= child_flags) &&
    any(parent_flags < child_flags)
}


# ============================================================
# 29. 可选：单个个体模拟测试
# ============================================================
# 这里不再默认使用 init_fit[1, ]。
# 改成使用 L0 最大的个体，更容易观察 rebound / reservoir driven dynamics。

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
  
  cat("\nSingle-subject CD8 killing rate check:\n")
  print(
    sim_test %>%
      summarise(
        id = first(id),
        group = first(group),
        CD8_count_G = first(CD8_count_G),
        k_CD8_G = first(k_CD8_G),
        CD8_killing_rate_G = first(CD8_killing_rate_G),
        L0 = first(L),
        I0 = first(I),
        V0 = first(V)
      )
  )
  
  p_sim_test <- sim_test %>%
    pivot_longer(
      cols = c(T, L, I, V, CD4_total),
      names_to = "state",
      values_to = "value"
    ) %>%
    ggplot(aes(time, log10(pmax(value, EPS)), color = state)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(
      title = paste0(
        "Single-subject simulation test: id = ",
        first(sim_test$id),
        ", CD8 killing rate = ",
        signif(first(sim_test$CD8_killing_rate_G), 3),
        " day^-1"
      ),
      x = "Time",
      y = "log10(value)"
    )
  
  print(p_sim_test)
  
  save_plot_png_pdf(
    p_sim_test,
    "Fig0_single_subject_simulation_test",
    width = 7,
    height = 5
  )
}


# ============================================================
# 30. 可选：先诊断 M0 初始参数下的残差
# ============================================================

flags_M0 <- c(
  za = 0,
  zp = 0,
  zk = 0
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
# 31. 可选：先测试 M0
# ============================================================

if (RUN_TEST_M0) {
  
  test_fit_M0 <- fit_one_model(
    model_row = candidate_models[1, ],
    start_free_pars = default_start_free_pars,
    start_id = 0
  )
  
  print(test_fit_M0$summary)
  print(test_fit_M0$params)
}


# ============================================================
# 32. 对 8 个候选模型执行 anchored random multi-start fitting
# ============================================================

random_fit_list <- vector("list", nrow(candidate_models))

for (i in seq_len(nrow(candidate_models))) {
  
  model_row_i <- candidate_models[i, ]
  
  anchor_starts_i <- list()
  
  if (i > 1) {
    
    for (j in seq_len(i - 1)) {
      
      parent_row_j <- candidate_models[j, ]
      
      if (is_nested_parent(parent_row_j, model_row_i)) {
        
        parent_params_j <- random_fit_list[[j]]$best_fit$params
        
        if (
          nrow(parent_params_j) > 0 &&
          !is.na(parent_params_j$a_L[1]) &&
          !is.na(parent_params_j$p[1]) &&
          !is.na(parent_params_j$k_CD8[1])
        ) {
          
          anchor_start_j <- params_row_to_start_free(parent_params_j[1, ])
          
          if (
            all(is.finite(anchor_start_j)) &&
            all(anchor_start_j[c("a_L", "p", "k_CD8")] > 0)
          ) {
            anchor_starts_i[[length(anchor_starts_i) + 1]] <- anchor_start_j
          }
        }
      }
    }
  }
  
  random_fit_list[[i]] <- fit_one_model_random_multistart(
    model_row = model_row_i,
    n_start = N_RANDOM_START,
    seed = RANDOM_SEED_BASE + i,
    anchor_starts = anchor_starts_i
  )
}


# ============================================================
# 33. 提取每个模型的 best-fit 结果
# ============================================================

fit_list_best <- purrr::map(
  random_fit_list,
  "best_fit"
)

model_fit_table_ms <- purrr::map_dfr(
  fit_list_best,
  "summary"
)

best_fit_params_ms <- purrr::map_dfr(
  fit_list_best,
  "params"
)


# ============================================================
# 34. 保存所有随机初值拟合记录
# ============================================================

all_random_fit_summaries <- purrr::map_dfr(
  seq_along(random_fit_list),
  function(i) {
    random_fit_list[[i]]$all_summaries
  }
)


# ============================================================
# 35. 基于 best-fit 计算 AICc 排名
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
# 36. nested RSS check
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
        parent_zp = zp,
        parent_zk = zk
      ),
    by = "parent"
  ) %>%
  left_join(
    candidate_models %>%
      rename(
        child = model,
        child_za = za,
        child_zp = zp,
        child_zk = zk
      ),
    by = "child"
  ) %>%
  filter(
    parent != child,
    parent_za <= child_za,
    parent_zp <= child_zp,
    parent_zk <= child_zk,
    parent_za < child_za |
      parent_zp < child_zp |
      parent_zk < child_zk
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
# 37. ΔAICc < 5 的 most parsimonious model set
# ============================================================

most_parsimonious_models <- aicc_table_ms %>%
  filter(delta_AICc < 5) %>%
  arrange(AICc)


# ============================================================
# 38. 参数边界检查
# ============================================================

boundary_check <- best_fit_params_ms %>%
  mutate(
    a_L_near_lower = !is.na(a_L) & a_L <= 1.01e-8,
    a_L_near_upper = !is.na(a_L) & a_L >= 0.99,
    p_near_lower = !is.na(p) & p <= 1.01,
    p_near_upper = !is.na(p) & p >= 0.99e8,
    k_CD8_near_lower = !is.na(k_CD8) & k_CD8 <= 1.01e-12,
    k_CD8_near_upper = !is.na(k_CD8) & k_CD8 >= 0.99e-5,
    theta_a_near_boundary = !is.na(theta_a) & abs(theta_a) >= 4.95,
    theta_p_near_boundary = !is.na(theta_p) & abs(theta_p) >= 4.95,
    theta_k_near_boundary = !is.na(theta_k) & abs(theta_k) >= 4.95
  )

print(boundary_check)

if (
  any(
    boundary_check$a_L_near_lower |
    boundary_check$a_L_near_upper |
    boundary_check$p_near_lower |
    boundary_check$p_near_upper |
    boundary_check$k_CD8_near_lower |
    boundary_check$k_CD8_near_upper |
    boundary_check$theta_a_near_boundary |
    boundary_check$theta_p_near_boundary |
    boundary_check$theta_k_near_boundary,
    na.rm = TRUE
  )
) {
  warning(
    "Some fitted parameters are close to bounds. ",
    "Interpret these models with caution."
  )
}


# ============================================================
# 39. 打印结果
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
# 40. 保存模型选择结果
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
  all_random_fit_summaries,
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
# 41. 保存最佳模型的残差表和预测表
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
    zp = as.numeric(best_model_row$zp),
    zk = as.numeric(best_model_row$zk)
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
}


# ============================================================
# 42. 保存本次运行设置
# ============================================================

run_settings <- tibble(
  item = c(
    "N_RANDOM_START",
    "RANDOM_SEED_BASE",
    "USE_VARIABLE_SCALING",
    "FIT_VARIABLES",
    "fixed_d_T",
    "fixed_beta",
    "fixed_f_L",
    "fixed_d_L",
    "fixed_delta0",
    "fixed_c",
    "k_CD8_start",
    "k_CD8_lower",
    "k_CD8_upper",
    "CD4_unit_assumption",
    "CD8_unit_assumption",
    "pVL_LOD_handling",
    "DNA_handling",
    "out_dir"
  ),
  value = c(
    as.character(N_RANDOM_START),
    as.character(RANDOM_SEED_BASE),
    as.character(USE_VARIABLE_SCALING),
    paste(FIT_VARIABLES, collapse = ","),
    as.character(fixed_pars["d_T"]),
    as.character(fixed_pars["beta"]),
    as.character(fixed_pars["f_L"]),
    as.character(fixed_pars["d_L"]),
    as.character(fixed_pars["delta0"]),
    as.character(fixed_pars["c"]),
    as.character(default_start_free_pars["k_CD8"]),
    "1e-12",
    "1e-5",
    "input CD4 is already cells/mL",
    "input CD8 is already cells/mL",
    "pVL value and LOD are raw copies/mL, both converted to log10; censoring judged on raw scale",
    "DNA excluded from fitting, residuals, RSS, AICc, observed-vs-predicted plots, and timecourse plots",
    out_dir
  )
)

write.csv(
  run_settings,
  file.path(out_dir, "run_settings.csv"),
  row.names = FALSE
)


# ============================================================
# 43. 可视化 1：AICc 模型排序
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
    title = "Model selection by AICc",
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
  width = 7,
  height = 5
)


# ============================================================
# 44. 可视化 2：最佳模型 observed vs predicted
# ============================================================

if (exists("best_residual_table") && !is.null(best_residual_table)) {
  
  p_obs_pred <- best_residual_table %>%
    ggplot(aes(x = value, y = pred)) +
    geom_point(aes(shape = is_censored), alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ variable, scales = "free") +
    theme_bw() +
    labs(
      title = paste0("Observed vs predicted values: ", best_model_name),
      x = "Observed value",
      y = "Predicted value",
      shape = "Censored"
    )
  
  print(p_obs_pred)
  
  save_plot_png_pdf(
    p_obs_pred,
    "Fig2_observed_vs_predicted_best_model",
    width = 9,
    height = 4.5
  )
}


# ============================================================
# 45. 可视化 3：最佳模型残差分布
# ============================================================

if (exists("best_residual_table") && !is.null(best_residual_table)) {
  
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
}


# ============================================================
# 46. 可视化 4：最佳模型时间序列拟合
# ============================================================

if (exists("best_residual_table") && !is.null(best_residual_table)) {
  
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
    height = 7
  )
}


# ============================================================
# 47. 可视化 5：最佳模型 early vs late 参数比较
# ============================================================

if (exists("best_fit_overall") && !is.null(best_fit_overall$params)) {
  
  best_param_plot_data <- best_fit_overall$params %>%
    select(
      model,
      a_L_late,
      a_L_early,
      p_late,
      p_early,
      k_CD8_late,
      k_CD8_early
    ) %>%
    pivot_longer(
      cols = -model,
      names_to = "parameter_group",
      values_to = "value"
    ) %>%
    mutate(
      parameter = case_when(
        grepl("^a_L", parameter_group) ~ "a_L",
        grepl("^p_", parameter_group) ~ "p",
        grepl("^k_CD8", parameter_group) ~ "k_CD8",
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
    "Fig5_best_model_parameter_comparison",
    width = 8,
    height = 4.5
  )
}


# ============================================================
# 48. 可视化 6：最佳模型不同变量的 RSS 贡献
# ============================================================

if (exists("best_residual_table") && !is.null(best_residual_table)) {
  
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
    "Fig6_weighted_RSS_by_variable",
    width = 7,
    height = 5
  )
}


# ============================================================
# 49. 可视化 7：最佳模型 CD8 effective killing rate
# ============================================================

if (exists("best_sim") && !is.null(best_sim)) {
  
  cd8_killing_summary <- best_sim %>%
    distinct(id, group, CD8_count_G, k_CD8_G, CD8_killing_rate_G) %>%
    arrange(group, id)
  
  write.csv(
    cd8_killing_summary,
    file.path(out_dir, "best_model_CD8_killing_rate_by_id.csv"),
    row.names = FALSE
  )
  
  p_cd8_kill <- cd8_killing_summary %>%
    ggplot(aes(x = group, y = CD8_killing_rate_G)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    theme_bw() +
    labs(
      title = paste0("Effective CD8 killing rate: ", best_model_name),
      x = "Group",
      y = expression(k[CD8] * " × CD8 count (day"^-1 * ")")
    )
  
  print(p_cd8_kill)
  
  save_plot_png_pdf(
    p_cd8_kill,
    "Fig7_effective_CD8_killing_rate",
    width = 7,
    height = 5
  )
}


# ============================================================
# 50. 最终提示
# ============================================================

cat("\nAll results saved to:\n")
cat(out_dir, "\n")

cat("\nFigures saved to:\n")
cat(fig_dir, "\n")

cat("\nMain files to check:\n")
cat("1. data_scale_check.csv\n")
cat("2. CD4_CD8_unit_check.csv\n")
cat("3. aicc_table_multistart.csv\n")
cat("4. best_fit_params_multistart.csv\n")
cat("5. nested_rss_check.csv\n")
cat("6. parameter_boundary_check.csv\n")
cat("7. best_model_residual_table.csv\n")
cat("8. best_model_RSS_by_variable.csv\n")
cat("9. best_model_CD8_killing_rate_by_id.csv\n")

cat("\nMain figures to show:\n")
cat("1. Fig1_AICc_model_ranking.png\n")
cat("2. Fig2_observed_vs_predicted_best_model.png\n")
cat("3. Fig3_residual_by_variable_best_model.png\n")
cat("4. Fig4_timecourse_observed_predicted_best_model.png\n")
cat("5. Fig5_best_model_parameter_comparison.png\n")
cat("6. Fig6_weighted_RSS_by_variable.png\n")
cat("7. Fig7_effective_CD8_killing_rate.png\n")

cat("\nBest model:\n")
cat(best_model_name, "\n")