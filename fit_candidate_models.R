library(deSolve)
library(tidyverse)
library(minpack.lm)

source("model_without_CD8.R")

data_clean <- read.csv("../Rdata/raw_data.csv")
init_table <- read.csv("../Rdata/init_table_I0_replaced.csv")


# ============================================================
# 1. 固定参数
# ============================================================

fixed_pars <- c(
  d_T = 0.01,
  beta = 1.58e-8,
  f_L = 0.005,
  d_L = 5e-4,
  c = 23
)


# ============================================================
# 2. 测试用参数
# ============================================================

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


# ============================================================
# 3. 初始拟合参数
# ============================================================

start_free_pars <- c(
  a_L = 1e-3,
  delta = 1.4,
  p = 6e3,
  theta_a = 0,
  theta_delta = 0,
  theta_p = 0
)


# ============================================================
# 4. 数据预处理
# ============================================================

valid_ids <- init_table %>%
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
  ) %>%
  pull(id)


data_fit <- data_clean %>%
  mutate(
    id = as.character(id),
    group = as.character(group),
    time = as.numeric(time),
    variable = as.character(variable),
    value = as.numeric(value),
    lod = as.numeric(lod)
  ) %>%
  filter(
    id %in% valid_ids,
    variable %in% c("pVL", "DNA", "CD4"),
    !is.na(time),
    !is.na(value)
  ) %>%
  mutate(
    value_raw = value,
    lod_raw = lod,
    value = case_when(
      variable == "pVL" ~ log10(pmax(value, 1e-12)),
      variable == "CD4" ~ log10(pmax(value, 1e-12)),
      variable == "DNA" ~ value,
      TRUE ~ value
    ),
    lod = case_when(
      variable == "pVL" & !is.na(lod) ~ log10(pmax(lod, 1e-12)),
      TRUE ~ lod
    )
  )


init_fit <- init_table %>%
  mutate(
    id = as.character(id),
    group = as.character(group),
    T0 = as.numeric(T0),
    L0 = as.numeric(L0),
    V0 = as.numeric(V0)
  ) %>%
  filter(id %in% valid_ids)


# ============================================================
# 5. 检查数据尺度
# ============================================================

scale_check <- data_fit %>%
  group_by(variable) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(scale_check)


# ============================================================
# 6. 8 个候选模型
# ============================================================

candidate_models <- tibble::tribble(
  ~model,        ~za, ~zdelta, ~zp,
  "M0",            0,       0,   0,
  "Ma",            1,       0,   0,
  "Mdelta",        0,       1,   0,
  "Mp",            0,       0,   1,
  "Ma_delta",      1,       1,   0,
  "Ma_p",          1,       0,   1,
  "Mdelta_p",      0,       1,   1,
  "Ma_delta_p",    1,       1,   1
)


# ============================================================
# 7. 单个个体模拟函数
# ============================================================

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


# ============================================================
# 8. 根据候选模型生成自由参数初值
# ============================================================

make_start_vector <- function(model_flags) {
  
  x0 <- c(
    log_a_L = log(unname(start_free_pars["a_L"])),
    log_delta = log(unname(start_free_pars["delta"])),
    log_p = log(unname(start_free_pars["p"]))
  )
  
  if (model_flags["za"] == 1) {
    x0 <- c(x0, theta_a = unname(start_free_pars["theta_a"]))
  }
  
  if (model_flags["zdelta"] == 1) {
    x0 <- c(x0, theta_delta = unname(start_free_pars["theta_delta"]))
  }
  
  if (model_flags["zp"] == 1) {
    x0 <- c(x0, theta_p = unname(start_free_pars["theta_p"]))
  }
  
  x0
}


# ============================================================
# 9. 参数下界
# ============================================================

make_lower_vector <- function(model_flags) {
  
  lower <- c(
    log_a_L = log(1e-8),
    log_delta = log(1e-3),
    log_p = log(1)
  )
  
  if (model_flags["za"] == 1) {
    lower <- c(lower, theta_a = -5)
  }
  
  if (model_flags["zdelta"] == 1) {
    lower <- c(lower, theta_delta = -5)
  }
  
  if (model_flags["zp"] == 1) {
    lower <- c(lower, theta_p = -5)
  }
  
  lower
}


# ============================================================
# 10. 参数上界
# ============================================================

make_upper_vector <- function(model_flags) {
  
  upper <- c(
    log_a_L = log(1),
    log_delta = log(10),
    log_p = log(1e7)
  )
  
  if (model_flags["za"] == 1) {
    upper <- c(upper, theta_a = 5)
  }
  
  if (model_flags["zdelta"] == 1) {
    upper <- c(upper, theta_delta = 5)
  }
  
  if (model_flags["zp"] == 1) {
    upper <- c(upper, theta_p = 5)
  }
  
  upper
}


# ============================================================
# 11. 把 nls.lm 的自由参数转回模型参数
# ============================================================

unpack_pars <- function(x, model_flags) {
  
  pars <- c(
    fixed_pars,
    a_L = exp(unname(x["log_a_L"])),
    delta = exp(unname(x["log_delta"])),
    p = exp(unname(x["log_p"])),
    theta_a = 0,
    theta_delta = 0,
    theta_p = 0
  )
  
  if (model_flags["za"] == 1) {
    pars["theta_a"] <- unname(x["theta_a"])
  }
  
  if (model_flags["zdelta"] == 1) {
    pars["theta_delta"] <- unname(x["theta_delta"])
  }
  
  if (model_flags["zp"] == 1) {
    pars["theta_p"] <- unname(x["theta_p"])
  }
  
  pars
}


# ============================================================
# 12. 模拟全部个体
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
# 13. 预测函数
# ============================================================

make_predictions <- function(sim) {
  
  sim %>%
    select(
      id,
      group,
      time,
      pred_pVL_log10,
      pred_DNA_log10,
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
        pred_DNA_log10 = "DNA",
        pred_CD4_log10 = "CD4"
      )
    )
}


# ============================================================
# 14. 残差函数
# ============================================================

make_residuals <- function(x, model_flags) {
  
  pars <- unpack_pars(x, model_flags)
  
  sim <- simulate_all(
    pars = pars,
    model_flags = model_flags
  )
  
  if (is.null(sim)) {
    return(rep(1e6, nrow(data_fit)))
  }
  
  pred <- make_predictions(sim)
  
  joined <- data_fit %>%
    left_join(
      pred,
      by = c("id", "group", "time", "variable")
    )
  
  if (any(is.na(joined$pred)) || any(!is.finite(joined$pred))) {
    return(rep(1e6, nrow(data_fit)))
  }
  
  resid <- joined$pred - joined$value
  
  if (any(!is.finite(resid))) {
    return(rep(1e6, nrow(data_fit)))
  }
  
  resid
}


# ============================================================
# 15. 诊断残差函数
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
        any_nonfinite = any(
          !is.finite(T) |
            !is.finite(L) |
            !is.finite(I) |
            !is.finite(V)
        )
      )
  )
  
  pred <- make_predictions(sim)
  
  cat("\nPrediction rows:", nrow(pred), "\n")
  
  joined <- data_fit %>%
    left_join(
      pred,
      by = c("id", "group", "time", "variable")
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
  
  if (any(!is.finite(joined$pred))) {
    cat("\nRows with non-finite predictions:\n")
    print(
      joined %>%
        filter(!is.finite(pred)) %>%
        select(id, group, time, variable, value, pred) %>%
        head(30)
    )
  }
  
  resid <- joined$pred - joined$value
  
  cat("\nResidual check:\n")
  cat("RSS:", sum(resid^2, na.rm = TRUE), "\n")
  cat("Any non-finite residual:", any(!is.finite(resid)), "\n")
  
  invisible(joined)
}


# ============================================================
# 16. 拟合单个候选模型
# ============================================================

fit_one_model <- function(model_row) {
  
  model_flags <- c(
    za = as.numeric(model_row$za),
    zdelta = as.numeric(model_row$zdelta),
    zp = as.numeric(model_row$zp)
  )
  
  x0 <- make_start_vector(model_flags)
  
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
      zdelta = model_row$zdelta,
      zp = model_row$zp,
      n = nrow(data_fit),
      q = length(x0),
      K = length(x0) + 1,
      RSS = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      convergence_info = NA_integer_
    )
    
    param_row <- tibble(
      model = model_row$model,
      a_L = NA_real_,
      delta = NA_real_,
      p = NA_real_,
      theta_a = NA_real_,
      theta_delta = NA_real_,
      theta_p = NA_real_,
      a_L_late = NA_real_,
      a_L_early = NA_real_,
      delta_late = NA_real_,
      delta_early = NA_real_,
      p_late = NA_real_,
      p_early = NA_real_
    )
    
    return(list(
      fit = NULL,
      summary = summary_row,
      params = param_row
    ))
  }
  
  x_hat <- fit$par
  pars_hat <- unpack_pars(x_hat, model_flags)
  
  resid_hat <- make_residuals(x_hat, model_flags)
  
  RSS <- sum(resid_hat^2)
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
  
  delta_late <- unname(pars_hat["delta"])
  delta_early <- unname(pars_hat["delta"]) *
    exp(unname(model_flags["zdelta"]) * unname(pars_hat["theta_delta"]))
  
  p_late <- unname(pars_hat["p"])
  p_early <- unname(pars_hat["p"]) *
    exp(-unname(model_flags["zp"]) * unname(pars_hat["theta_p"]))
  
  summary_row <- tibble(
    model = model_row$model,
    za = model_row$za,
    zdelta = model_row$zdelta,
    zp = model_row$zp,
    n = n,
    q = q,
    K = K,
    RSS = RSS,
    logLik = logLik,
    AIC = AIC,
    AICc = AICc,
    convergence_info = fit$info
  )
  
  param_row <- tibble(
    model = model_row$model,
    a_L = unname(pars_hat["a_L"]),
    delta = unname(pars_hat["delta"]),
    p = unname(pars_hat["p"]),
    theta_a = unname(pars_hat["theta_a"]),
    theta_delta = unname(pars_hat["theta_delta"]),
    theta_p = unname(pars_hat["theta_p"]),
    a_L_late = a_L_late,
    a_L_early = a_L_early,
    delta_late = delta_late,
    delta_early = delta_early,
    p_late = p_late,
    p_early = p_early
  )
  
  list(
    fit = fit,
    summary = summary_row,
    params = param_row
  )
}


# ============================================================
# 17. 随机初值生成函数
# ============================================================

make_random_start <- function() {
  
  c(
    a_L = exp(runif(1, log(1e-5), log(1e-1))),
    delta = exp(runif(1, log(0.1), log(5))),
    p = exp(runif(1, log(1e2), log(1e5))),
    theta_a = runif(1, -2, 2),
    theta_delta = runif(1, -2, 2),
    theta_p = runif(1, -2, 2)
  )
}


make_random_start_set <- function(n_start = 10, seed = 123) {
  
  set.seed(seed)
  
  start_set <- vector("list", n_start)
  
  for (i in seq_len(n_start)) {
    start_set[[i]] <- make_random_start()
  }
  
  start_set
}


# ============================================================
# 18. 指定初值拟合单个模型
# ============================================================

fit_one_model_from_start <- function(model_row, start_pars_i, start_id) {
  
  old_start_free_pars <- start_free_pars
  
  on.exit({
    start_free_pars <<- old_start_free_pars
  }, add = TRUE)
  
  start_free_pars <<- start_pars_i
  
  fit_i <- fit_one_model(model_row)
  
  fit_i$summary <- fit_i$summary %>%
    mutate(start_id = start_id)
  
  fit_i$params <- fit_i$params %>%
    mutate(start_id = start_id)
  
  fit_i
}


# ============================================================
# 19. 单个候选模型的 multi-start fitting
# ============================================================

fit_one_model_random_multistart <- function(model_row, n_start = 10, seed = 123) {
  
  start_set <- make_random_start_set(
    n_start = n_start,
    seed = seed
  )
  
  fit_candidates <- purrr::map(
    seq_along(start_set),
    function(i) {
      cat(
        "\nModel:", model_row$model,
        "| random start:", i,
        "\n"
      )
      
      fit_one_model_from_start(
        model_row = model_row,
        start_pars_i = start_set[[i]],
        start_id = i
      )
    }
  )
  
  summary_candidates <- purrr::map_dfr(
    fit_candidates,
    "summary"
  )
  
  valid_candidates <- summary_candidates %>%
    filter(!is.na(logLik), is.finite(logLik))
  
  if (nrow(valid_candidates) == 0) {
    
    failed_summary <- summary_candidates %>%
      slice(1) %>%
      mutate(
        best_start_id = NA_integer_,
        n_random_start = n_start
      )
    
    failed_params <- purrr::map_dfr(
      fit_candidates,
      "params"
    ) %>%
      slice(1) %>%
      mutate(
        best_start_id = NA_integer_,
        n_random_start = n_start
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
      n_random_start = n_start
    )
  
  best_fit$params <- best_fit$params %>%
    mutate(
      best_start_id = best_start_id,
      n_random_start = n_start
    )
  
  list(
    best_fit = best_fit,
    all_summaries = summary_candidates
  )
}


# ============================================================
# 20. 可选：单个个体模拟测试
# ============================================================

times <- seq(0, max(data_clean$time, na.rm = TRUE), by = 0.1)

sim_test <- simulate_one(
  id_row = init_table[1, ],
  pars = test_pars,
  times = times,
  model_flags = test_flags
)

if (is.null(sim_test)) {
  stop("Simulation failed.")
}

print(
  sim_test %>%
    pivot_longer(
      cols = c(T, L, I, V, CD4_total),
      names_to = "state",
      values_to = "value"
    ) %>%
    ggplot(aes(time, log10(pmax(value, 1e-12)), color = state)) +
    geom_line() +
    theme_bw()
)


# ============================================================
# 21. 可选：先诊断 M0 初始参数下的残差
# ============================================================

flags_M0 <- c(
  za = 0,
  zdelta = 0,
  zp = 0
)

x0_M0 <- make_start_vector(flags_M0)

debug_joined_M0 <- debug_residuals(
  x = x0_M0,
  model_flags = flags_M0
)


# ============================================================
# 22. 可选：先测试 M0
# ============================================================

test_fit_M0 <- fit_one_model(candidate_models[1, ])

print(test_fit_M0$summary)
print(test_fit_M0$params)


# ============================================================
# 23. 对 8 个候选模型执行 random multi-start fitting
# ============================================================

random_fit_list <- purrr::map(
  seq_len(nrow(candidate_models)),
  function(i) {
    fit_one_model_random_multistart(
      model_row = candidate_models[i, ],
      n_start = 5,
      seed = 1000 + i
    )
  }
)


# ============================================================
# 24. 提取每个模型的 best-fit 结果
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
# 25. 保存所有随机初值拟合记录
# ============================================================

all_random_fit_summaries <- purrr::map_dfr(
  seq_along(random_fit_list),
  function(i) {
    random_fit_list[[i]]$all_summaries
  }
)


# ============================================================
# 26. 基于 random multi-start best-fit 计算 AICc 排名
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
# 27. ΔAICc < 5 的 most parsimonious model set
# ============================================================

most_parsimonious_models <- aicc_table_ms %>%
  filter(delta_AICc < 5) %>%
  arrange(AICc)


# ============================================================
# 28. 打印结果
# ============================================================

print(best_fit_params_ms)

aicc_table_ms %>%
  select(
    rank,
    model,
    RSS,
    logLik,
    AICc,
    delta_AICc,
    Akaike_weight,
    q,
    K,
    best_start_id,
    n_random_start,
    convergence_info
  ) %>%
  print(n = Inf)

most_parsimonious_models %>%
  select(
    rank,
    model,
    RSS,
    logLik,
    AICc,
    delta_AICc,
    Akaike_weight,
    q,
    K
  ) %>%
  print(n = Inf)


# ============================================================
# 29. 保存 random multi-start 模型选择结果
# ============================================================

dir.create(
  "../Rdata/model_selection_results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  best_fit_params_ms,
  "../Rdata/model_selection_results/best_fit_params_multistart.csv",
  row.names = FALSE
)

write.csv(
  model_fit_table_ms,
  "../Rdata/model_selection_results/model_fit_summary_multistart.csv",
  row.names = FALSE
)

write.csv(
  aicc_table_ms,
  "../Rdata/model_selection_results/aicc_table_multistart.csv",
  row.names = FALSE
)

write.csv(
  most_parsimonious_models,
  "../Rdata/model_selection_results/most_parsimonious_models_deltaAICc_lt5.csv",
  row.names = FALSE
)

write.csv(
  all_random_fit_summaries,
  "../Rdata/model_selection_results/all_random_fit_summaries.csv",
  row.names = FALSE
)