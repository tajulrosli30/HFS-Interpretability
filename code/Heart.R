suppressPackageStartupMessages({
  library(FuzzyR)
  library(tidyverse)
  library(readr)
  library(janitor)
  library(caret)
})

# -------------------------------
# 0) CONFIG & DATA
# -------------------------------
data_path <- "G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/Data/heart_disease_uci.csv"

# Correct target column name from CSV
target_col <- "num"
feat_cols  <- c("age", "trestbps", "chol")
set.seed(42)

# Load and clean data
df <- read_csv(data_path, show_col_types = FALSE) |>
  clean_names() |>
  select(all_of(c(feat_cols, target_col))) |>
  drop_na()

# Convert target to binary (0 / 1)
df <- df %>% mutate(!!target_col := ifelse(.data[[target_col]] %in% c(1, 2, 3, 4, "1", "yes", "Yes"), 1L, 0L))

# Train/test split
trainIndex <- createDataPartition(df[[target_col]], p = 0.8, list = FALSE)
train <- df[trainIndex, ]
test  <- df[-trainIndex, ]

# Ensure numeric features
train <- train %>% mutate(across(all_of(feat_cols), as.numeric))
test  <- test  %>% mutate(across(all_of(feat_cols), as.numeric))

# -------------------------------
# Helper functions
# -------------------------------
gauss_eval <- function(x, sigma, c) exp(-0.5 * ((x - c)/sigma)^2)

mf_argmax <- function(x, mf_list) {
  degs <- sapply(mf_list, function(m) gauss_eval(x, m$params[1], m$params[2]))
  which.max(degs)
}

build_feature_mfs <- function(x, k=3) {
  r <- range(x, na.rm = TRUE)
  centers <- seq(r[1], r[2], length.out = k)
  sigmas <- rep((r[2]-r[1])/(2*k), k)
  mfs <- vector("list", k)
  for (i in seq_len(k)) mfs[[i]] <- list(type = "gaussmf", params = c(sigmas[i], centers[i]))
  mfs
}

safe_eval <- function(x, fis) {
  res <- evalfis(as.matrix(x), fis)
  res[is.na(res)] <- 0.5
  pmax(0, pmin(1, res))
}

# -------------------------------
# 1) FLS – Full rule base (3-in)
# -------------------------------
build_fls_full <- function(feat_cols, train, target_col) {
  fis <- newfis("FLS_full", fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  mf_defs <- list()
  for (i in seq_along(feat_cols)) {
    x <- train[[feat_cols[i]]]
    fis <- addvar(fis, "input", feat_cols[i], range(x, na.rm = TRUE))
    mfs <- build_feature_mfs(x, 3)
    mf_defs[[i]] <- mfs
    for (j in seq_along(mfs)) {
      fis <- addmf(fis, "input", i, paste0("G", j), "gaussmf", mfs[[j]]$params)
    }
  }
  fis <- addvar(fis, "output", "Risk", c(0,1))
  fis <- addmf(fis, "output", 1, "NoDisease", "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Disease",   "trimf", c(0.5,1,1))
  
  grid <- expand.grid(rep(list(1:3), length(feat_cols)))
  rule_mat <- as.matrix(cbind(grid, 1, 1, 1))
  fis <- addrule(fis, rule_mat)
  
  list(fis = fis, n_rules_total = nrow(rule_mat))
}

fls <- build_fls_full(feat_cols, train, target_col)
preds_fls <- safe_eval(test[, feat_cols], fls$fis)
yhat_fls <- ifelse(preds_fls > 0.5, 1, 0)
res_fls <- list(conf = table(Predicted = yhat_fls, Actual = test[[target_col]]),
                acc = mean(yhat_fls == test[[target_col]]))

cat("\n=== FLS (Full Rule Base) ===\n")
print(res_fls$conf)
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Total Rules:", fls$n_rules_total, "\n")

# -------------------------------
# 2) HFS Serial (2 subsystems, 2 layers) – Data-driven
# -------------------------------
build_2iso_sparse <- function(df, in1, in2, target_col) {
  fis <- newfis(paste0("FIS_", in1, "_", in2), fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  feat_pair <- c(in1, in2)
  mf_defs <- list()
  for (var_idx in seq_along(feat_pair)) {
    vec <- df[[feat_pair[var_idx]]]
    fis <- addvar(fis, "input", feat_pair[var_idx], range(vec, na.rm = TRUE))
    mfs <- build_feature_mfs(vec, 3)
    mf_defs[[var_idx]] <- mfs
    for (j in seq_along(mfs)) {
      fis <- addmf(fis, "input", var_idx, paste0("G", j), "gaussmf", mfs[[j]]$params)
    }
  }
  fis <- addvar(fis, "output", "Risk", c(0,1))
  fis <- addmf(fis, "output", 1, "NoDisease", "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Disease",   "trimf", c(0.5,1,1))
  
  idx_mat <- t(apply(as.matrix(df[, feat_pair]), 1, function(r) {
    sapply(seq_along(feat_pair), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  idx_df <- as.data.frame(idx_mat)
  colnames(idx_df) <- c("V1","V2")
  
  rule_tbl <- bind_cols(idx_df, tibble(cls = df[[target_col]])) %>%
    count(V1, V2, cls, name = "n") %>%
    group_by(V1, V2) %>%
    slice_max(order_by = n, n = 1) %>%
    ungroup() %>%
    mutate(out_mf = ifelse(cls == 0, 1L, 2L), weight = 1)
  
  rule_mat <- as.matrix(cbind(rule_tbl$V1, rule_tbl$V2,
                              rule_tbl$out_mf, rule_tbl$weight, 1L))
  fis <- addrule(fis, rule_mat)
  
  list(fis = fis, n_rules_total = nrow(rule_mat))
}

# S1: (age, trestbps) -> risk1
S1_ser <- build_2iso_sparse(train, "age", "trestbps", target_col)
train$risk1 <- safe_eval(train[, c("age","trestbps")], S1_ser$fis)
test$risk1  <- safe_eval(test[, c("age","trestbps")], S1_ser$fis)

# S2: (risk1, chol) -> final risk
S2_ser <- build_2iso_sparse(train, "risk1", "chol", target_col)
pred_ser <- safe_eval(test[, c("risk1","chol")], S2_ser$fis)
yhat_ser <- ifelse(pred_ser > 0.5, 1, 0)
conf_ser <- table(Predicted = yhat_ser, Actual = test[[target_col]])
acc_ser  <- mean(yhat_ser == test[[target_col]])

cat("\n=== HFS Serial (2 subsystems) ===\n")
print(conf_ser)
cat("Accuracy:", sprintf("%.4f", acc_ser), "\n")
rules_ser <- S1_ser$n_rules_total + S2_ser$n_rules_total
cat("Total Rules:", rules_ser, "\n")

#######################################################################

# --- helpers: get numeric indices from your existing functions ---

source("G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/R code/Nauck_Index.R")
source("G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/R code/Fuzzy_index_V2.R")

# ==== 1) Numeric wrappers around your printed indices ====

nauck_index_num <- function(fis) {
  out <- tryCatch(capture.output(Nauck_I(fis)), error = function(e) return(character()))
  line <- tail(grep("Nauck's Index", out, value = TRUE), 1)
  num  <- suppressWarnings(as.numeric(sub(".*Nauck's Index\\s*:\\s*", "", line)))
  if (is.na(num)) {
    nums <- unlist(regmatches(out, gregexpr("[-+]?[0-9]*\\.?[0-9]+", out)))
    num  <- suppressWarnings(as.numeric(tail(nums, 1)))
  }
  if (is.na(num)) NA_real_ else num
}

fuzzy_index_num <- function(fis) {
  out <- tryCatch(capture.output(Fuzzy_Index(fis)), error = function(e) return(character()))
  line <- tail(grep("Fuzzy Index", out, value = TRUE), 1)
  num  <- suppressWarnings(as.numeric(sub(".*Fuzzy Index\\s*:\\s*", "", line)))
  if (is.na(num)) {
    nums <- unlist(regmatches(out, gregexpr("[-+]?[0-9]*\\.?[0-9]+", out)))
    num  <- suppressWarnings(as.numeric(tail(nums, 1)))
  }
  if (is.na(num)) NA_real_ else num
}

# ==== 2) Sanitize FIS: fix NA/non-positive Gaussian sigmas ====

fix_gauss_sigmas <- function(fis) {
  for (i in seq_along(fis$input)) {
    rng   <- fis$input[[i]]$range
    width <- diff(rng)
    mfs   <- fis$input[[i]]$mf
    if (is.null(mfs)) next
    k <- length(mfs)
    for (j in seq_along(mfs)) {
      m <- mfs[[j]]
      if (!identical(m$type, "gaussmf")) next
      sig <- m$params[1]
      cen <- m$params[2]
      invalid <- is.na(sig) || !is.finite(sig) || sig <= 0
      if (invalid) {
        # neighbor-based fallback; then floor by global range/(5k)
        left  <- if (j > 1 && identical(mfs[[j-1]]$type, "gaussmf")) abs(cen - mfs[[j-1]]$params[2]) else NA_real_
        right <- if (j < k && identical(mfs[[j+1]]$type, "gaussmf")) abs(mfs[[j+1]]$params[2] - cen) else NA_real_
        cand <- suppressWarnings(stats::median(c(left, right), na.rm = TRUE))
        if (!is.finite(cand)) cand <- width / max(2 * k, 1)
        fis$input[[i]]$mf[[j]]$params[1] <- max(cand, width / max(5 * k, 1), 1e-6)
      }
    }
  }
  fis
}

fix_layers <- function(layers) {
  lapply(layers, function(subs) lapply(subs, fix_gauss_sigmas))
}

# ==== 3) E_{jk} and H_mean (equal layer weights, min within layer) ====

E_from_fis <- function(fis, w_nauck = 0.5) {
  n <- tryCatch(nauck_index_num(fis), error = function(e) NA_real_)
  f <- tryCatch(fuzzy_index_num(fis),  error = function(e) NA_real_)
  if (is.na(n) || is.na(f)) return(NA_real_)
  w_nauck * n + (1 - w_nauck) * f
}

H_mean_equal <- function(layers, w_nauck = 0.5) {
  # drop NULL/empty layers
  layers <- layers[ vapply(layers, function(x) length(x) > 0, logical(1)) ]
  if (!length(layers)) return(NA_real_)
  # per-layer score = min over its subsystems
  layer_vals <- vapply(layers, function(fis_list) {
    Es <- vapply(fis_list, function(f) E_from_fis(f, w_nauck), numeric(1))
    Es <- Es[is.finite(Es)]
    if (!length(Es)) return(NA_real_) else min(Es)
  }, numeric(1))
  layer_vals <- layer_vals[is.finite(layer_vals)]
  if (!length(layer_vals)) return(NA_real_)
  mean(layer_vals)  # equal weights l_j = 1/q
}

# ==== 4) Build your layer lists (as you already do) ====

# FLS: 1 layer, 1 subsystem
layers_fls <- list(
  list(fls$fis)
)

# HFS Serial: 2 layers, 1 subsystem each
layers_ser <- list(
  list(S1_ser$fis),   # Layer 1: (age, trestbps)
  list(S2_ser$fis)    # Layer 2: (risk1, chol)
)

# HFS Parallel: NOT AVAILABLE for Heart
# We can set it to NULL or an empty list if you don't want to compute for it
layers_par <- NULL

# ==== 5) Sanitize all layers to remove NA sigmas BEFORE indices ====
layers_fls <- fix_layers(layers_fls)
layers_ser <- fix_layers(layers_ser)
layers_par <- fix_layers(layers_par)

# (Optional) quick assert: ensure no NA gauss sigma remains
.check_na_sigma <- function(fis) {
  any(sapply(fis$input, function(inp)
    any(sapply(inp$mf, function(m)
      identical(m$type, "gaussmf") && (is.na(m$params[1]) || !is.finite(m$params[1]) || m$params[1] <= 0)))))
}
stopifnot(!any(unlist(lapply(layers_fls, function(L) sapply(L, .check_na_sigma)))))
stopifnot(!any(unlist(lapply(layers_ser, function(L) sapply(L, .check_na_sigma)))))
stopifnot(!any(unlist(lapply(layers_par, function(L) sapply(L, .check_na_sigma)))))

# ==== 6) Compute H_mean (equal Nauck/Fuzzy weight = 0.5) ====
H_fls <- H_mean_equal(layers_fls, w_nauck = 0.5)
H_ser <- H_mean_equal(layers_ser, w_nauck = 0.5)
H_par <- H_mean_equal(layers_par, w_nauck = 0.5)

cat(sprintf("\nH_mean — FLS: %.4f | HFS Serial: %.4f | HFS Parallel: %.4f\n", H_fls, H_ser, H_par))


