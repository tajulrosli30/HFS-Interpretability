#############################################
# CKD Dataset — FLS, HFS-Parallel, HFS-Serial
#############################################

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
data_path  <- "G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/Data/kidney_disease.csv"
target_col <- "classification"
feat_cols  <- c("age", "bp", "sg", "bgr", "sc")
set.seed(42)

df <- read_csv(data_path, show_col_types = FALSE) |>
  clean_names() |>
  select(all_of(c(feat_cols, target_col))) |>
  drop_na()

# Encode target as numeric 0/1
df <- df %>%
  mutate(!!target_col := ifelse(!!sym(target_col) == "ckd", 1L, 0L))

# Stratified train/test 80/20
trainIndex <- createDataPartition(df[[target_col]], p = 0.8, list = FALSE)
train <- df[trainIndex, ]
test  <- df[-trainIndex, ]

# Ensure numeric inputs
train <- train %>% mutate(across(all_of(feat_cols), as.numeric))
test  <- test %>% mutate(across(all_of(feat_cols), as.numeric))

# -------------------------------
# 1) Helpers
# -------------------------------
infer_nmf <- function(x) {
  ux <- sort(unique(x))
  if (length(ux) <= 3 || all(ux %in% c(0,1))) 2L else 3L
}
gauss_eval <- function(x, sigma, c) exp(-0.5 * ((x - c)/sigma)^2)
mf_argmax <- function(x, mf_list) {
  degs <- sapply(mf_list, function(m) gauss_eval(x, m$params[1], m$params[2]))
  which.max(degs)
}
crisp_to_class <- function(y) { ifelse(y < 0.5, 0, 1) }
km_centers_sigma <- function(x, k) {
  x <- as.numeric(x)
  if (sd(x, na.rm = TRUE) == 0 || length(unique(x)) < k) {
    c_min <- min(x, na.rm = TRUE); c_max <- max(x, na.rm = TRUE)
    centers <- if (k == 2) c(c_min, c_max) else seq(c_min, c_max, length.out = k)
  } else {
    km <- kmeans(x, centers = k, nstart = 20)
    centers <- sort(as.numeric(km$centers))
  }
  gaps <- diff(centers)
  sigmas <- if (length(gaps) > 0) pmax(gaps, diff(range(x, na.rm = TRUE))/(5*k)) else rep(diff(range(x))/(2*k), k)
  list(centers = centers, sigmas = sigmas)
}
build_feature_mfs <- function(x, k) {
  ini <- km_centers_sigma(x, k)
  centers <- ini$centers; sigmas <- ini$sigmas
  mfs <- vector("list", length = k)
  for (i in seq_len(k)) mfs[[i]] <- list(type = "gaussmf", params = c(sigmas[i], centers[i]))
  mfs
}
add_default_rule <- function(rule_mat, n_in) {
  r <- rep(1L, n_in)
  cbind(matrix(r, nrow = 1), output = 1L, weight = 0.1, connector = 1L)
}
safe_eval <- function(x, fis) {
  res <- evalfis(as.matrix(x), fis)
  res[is.na(res)] <- 0.5
  pmax(0, pmin(1, res))
}

# -------------------------------
# 2) Build FLS (5-in)
# -------------------------------
build_fls <- function(train, feat_cols, target_col, n_mf_vec) {
  fis <- newfis("FLS_model", fisType = "mamdani",
                andMethod = "min", orMethod = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  
  mf_defs <- vector("list", length(feat_cols))
  
  for (i in seq_along(feat_cols)) {
    x <- train[[feat_cols[i]]]
    fis <- addvar(fis, "input", feat_cols[i], range(x, na.rm = TRUE))
    mfs <- build_feature_mfs(x, n_mf_vec[i])
    mf_defs[[i]] <- mfs
    for (j in seq_along(mfs)) {
      fis <- addmf(fis, "input", i, paste0("G", j), "gaussmf", mfs[[j]]$params)
    }
  }
  
  fis <- addvar(fis, "output", "CKD_Risk", c(0, 1))
  fis <- addmf(fis, "output", 1, "NotCKD", "trimf", c(0, 0, 0.5))
  fis <- addmf(fis, "output", 1, "CKD",    "trimf", c(0.5, 1, 1))
  
  idx_mat <- t(apply(as.matrix(train[, feat_cols]), 1, function(r) {
    sapply(seq_along(feat_cols), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  
  idx_df <- as.data.frame(idx_mat)
  colnames(idx_df) <- paste0("V", seq_along(feat_cols))
  
  dat <- bind_cols(idx_df, tibble(cls = train[[target_col]]))
  
  rule_tbl <- dat %>%
    count(across(starts_with("V")), cls, name = "n") %>%
    group_by(across(starts_with("V"))) %>%
    slice_max(order_by = n, n = 1) %>%
    ungroup() %>%
    mutate(out_mf = ifelse(cls == 0, 1L, 2L),
           weight = pmax(0.5, n / max(n)))
  
  rule_mat <- as.matrix(cbind(as.matrix(rule_tbl[, paste0("V", seq_along(feat_cols))]),
                              rule_tbl$out_mf, rule_tbl$weight, 1L))
  
  fis <- addrule(fis, rbind(rule_mat, add_default_rule(rule_mat, length(feat_cols))))
  
  list(fis = fis, n_rules_total = nrow(rule_mat) + 1)
}

n_mf_vec <- sapply(train[feat_cols], infer_nmf)
fls <- build_fls(train, feat_cols, target_col, n_mf_vec)

eval_model <- function(fis, newdata, feat_cols, y_true) {
  preds <- safe_eval(newdata[, feat_cols], fis)
  y_hat <- ifelse(preds > 0.5, 1, 0)
  list(conf = table(Predicted = y_hat, Actual = y_true),
       acc = mean(y_hat == y_true))
}

res_fls <- eval_model(fls$fis, test, feat_cols, test[[target_col]])

cat("\n=== FLS (5-in) ===\n")
print(res_fls$conf)
cat("Accuracy:", res_fls$acc, "\n")
cat("Rules:", fls$n_rules_total, "\n")

# -------------------------------
# 3) Build 2-input FIS helper
# -------------------------------
build_2iso <- function(df, in1, in2, target_col) {
  fis <- newfis(paste0("FIS_", in1, "_", in2),
                fisType = "mamdani",
                andMethod = "min", orMethod = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  
  feat_pair <- c(in1, in2)
  n_mf_vec <- sapply(df[feat_pair], infer_nmf)
  
  mf_defs <- list()
  for (var_idx in seq_along(feat_pair)) {
    var_name <- feat_pair[var_idx]
    vec <- df[[var_name]]
    fis <- addvar(fis, "input", var_name, range(vec, na.rm = TRUE))
    mfs <- build_feature_mfs(vec, n_mf_vec[var_idx])
    mf_defs[[var_idx]] <- mfs
    for (j in seq_along(mfs)) {
      fis <- addmf(fis, "input", var_idx, paste0("G", j), "gaussmf", mfs[[j]]$params)
    }
  }
  
  fis <- addvar(fis, "output", "Risk", c(0, 1))
  fis <- addmf(fis, "output", 1, "NotCKD", "trimf", c(0, 0, 0.5))
  fis <- addmf(fis, "output", 1, "CKD",    "trimf", c(0.5, 1, 1))
  
  idx_mat <- t(apply(as.matrix(df[, feat_pair]), 1, function(r) {
    sapply(seq_along(feat_pair), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  
  idx_df <- as.data.frame(idx_mat)
  colnames(idx_df) <- paste0("V", seq_along(feat_pair))
  
  rule_tbl <- bind_cols(idx_df, tibble(cls = df[[target_col]])) %>%
    count(across(starts_with("V")), cls, name = "n") %>%
    group_by(across(starts_with("V"))) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(out_mf = ifelse(cls == 0, 1L, 2L),
           weight = pmax(0.5, n / max(n)))
  
  rule_mat <- as.matrix(cbind(rule_tbl[, paste0("V", seq_along(feat_pair))],
                              rule_tbl$out_mf, rule_tbl$weight, 1L))
  
  fis <- addrule(fis, rbind(rule_mat, add_default_rule(rule_mat, 2)))
  
  list(fis = fis, n_rules_total = nrow(rule_mat) + 1)
}

# -------------------------------
# 4) HFS Parallel (Updated Design)
# -------------------------------
# Layer 1
S1 <- build_2iso(train, "age", "bp", target_col)           # demographic_risk
S2 <- build_2iso(train, "sg", "bgr", target_col)           # metabolic_hydration_risk
train$y1 <- safe_eval(train[, c("age", "bp")], S1$fis)
train$y2 <- safe_eval(train[, c("sg", "bgr")], S2$fis)
test$y1  <- safe_eval(test[,  c("age", "bp")], S1$fis)
test$y2  <- safe_eval(test[,  c("sg", "bgr")], S2$fis)

# Layer 2
S3 <- build_2iso(train, "y1", "sc", target_col)            # renal_demographic_risk
train$y3 <- safe_eval(train[, c("y1", "sc")], S3$fis)
test$y3  <- safe_eval(test[,  c("y1", "sc")], S3$fis)

# Layer 3
S4 <- build_2iso(train, "y2", "y3", target_col)            # CKD Risk
pred_parallel <- safe_eval(test[, c("y2", "y3")], S4$fis)
pred_class <- ifelse(pred_parallel > 0.5, 1, 0)
acc_parallel <- mean(pred_class == test[[target_col]], na.rm = TRUE)
conf_parallel <- table(Predicted = pred_class, Actual = test[[target_col]])

# -------------------------------
# 5) HFS Serial
# -------------------------------
SS1 <- build_2iso(train, "age", "bp", target_col)
train$z1 <- safe_eval(train[, c("age", "bp")], SS1$fis)
test$z1  <- safe_eval(test[,  c("age", "bp")], SS1$fis)

SS2 <- build_2iso(train, "z1", "sg", target_col)
train$z2 <- safe_eval(train[, c("z1", "sg")], SS2$fis)
test$z2  <- safe_eval(test[,  c("z1", "sg")], SS2$fis)

SS3 <- build_2iso(train, "z2", "bgr", target_col)
train$z3 <- safe_eval(train[, c("z2", "bgr")], SS3$fis)
test$z3  <- safe_eval(test[,  c("z2", "bgr")], SS3$fis)

SS4 <- build_2iso(train, "z3", "sc", target_col)
pred_serial <- safe_eval(test[, c("z3", "sc")], SS4$fis)
pred_class_s <- ifelse(pred_serial > 0.5, 1, 0)
acc_serial <- mean(pred_class_s == test[[target_col]], na.rm = TRUE)
conf_serial <- table(Predicted = pred_class_s, Actual = test[[target_col]])

# -------------------------------
# 6) Results
# -------------------------------
cat("\n=== HFS Parallel ===\n")
print(conf_parallel)
cat("Accuracy:", acc_parallel, "\n")

cat("\n=== HFS Serial ===\n")
print(conf_serial)
cat("Accuracy:", acc_serial, "\n")

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

# ==== 4) Build your layer lists for Kidney ====

# FLS: 1 layer, 1 subsystem
layers_fls <- list( list(fls$fis) )

# HFS Serial: 4 layers, 1 subsystem each
layers_ser <- list(
  list(SS1$fis),
  list(SS2$fis),
  list(SS3$fis),
  list(SS4$fis)
)

# HFS Parallel: L1 has 2 subsystems, then 1 each for L2 & L3
layers_par <- list(
  list(S1$fis, S2$fis),  # Layer 1
  list(S3$fis),          # Layer 2
  list(S4$fis)           # Layer 3
)


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
