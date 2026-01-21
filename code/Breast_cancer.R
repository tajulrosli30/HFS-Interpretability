#############################################
# Breast Cancer — FLS, HFS Parallel, HFS Serial
# Binary target: diagnosis (B=0, M=1)
# Inputs (6): radius_mean, texture_mean, perimeter_mean,
#             area_mean, smoothness_mean, compactness_mean
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
# Change this to your file path if needed
data_path <- "G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/Data/data_breast_cancer.csv"
target_col <- "diagnosis"
feat_cols  <- c("radius_mean","texture_mean","perimeter_mean",
                "area_mean","smoothness_mean","compactness_mean")
set.seed(42)

df <- read_csv(data_path, show_col_types = FALSE) |>
  clean_names()

# harmonize column names to our expected ones if necessary
name_map <- c(
  "radius_mean"      = "radius_mean",
  "texture_mean"     = "texture_mean",
  "perimeter_mean"   = "perimeter_mean",
  "area_mean"        = "area_mean",
  "smoothness_mean"  = "smoothness_mean",
  "compactness_mean" = "compactness_mean",
  "diagnosis"        = "diagnosis"
)
# keep only needed cols (robust to different original ordering)
df <- df |> rename(any_of(name_map)) |> select(all_of(c(feat_cols, target_col))) |> drop_na()

# Encode target: B=0, M=1
df <- df %>% mutate(!!target_col := ifelse(.data[[target_col]] %in% c("M","malignant","Malignant",1, "1"), 1L, 0L))

# Stratified train/test 80/20
trainIndex <- createDataPartition(df[[target_col]], p = 0.8, list = FALSE)
train <- df[trainIndex, ]
test  <- df[-trainIndex, ]

# Ensure numeric inputs
train <- train %>% mutate(across(all_of(feat_cols), as.numeric))
test  <- test  %>% mutate(across(all_of(feat_cols), as.numeric))

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
to_latex_mat <- function(tab) {
  m <- as.matrix(tab)
  if (nrow(m) == 0 || ncol(m) == 0) return("\\begin{bmatrix}\\end{bmatrix}")
  rows <- apply(m, 1, function(r) paste0(r, collapse = " & "))
  paste0("\\begin{bmatrix}", paste(rows, collapse = " \\\\ "), " \\end{bmatrix}")
}

# -------------------------------
# 2) FLS (6-in) – data-driven rules
# -------------------------------
build_fls <- function(train, feat_cols, target_col, n_mf_vec) {
  fis <- newfis("FLS_6in", fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
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
  
  # Binary output
  fis <- addvar(fis, "output", "CancerRisk", c(0,1))
  fis <- addmf(fis, "output", 1, "Benign",     "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Malignant",  "trimf", c(0.5,1,1))
  
  # Map each row to MF indices
  idx_mat <- t(apply(as.matrix(train[, feat_cols]), 1, function(r) {
    sapply(seq_along(feat_cols), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  idx_df <- as.data.frame(idx_mat)
  colnames(idx_df) <- paste0("V", seq_along(feat_cols))
  
  dat <- bind_cols(idx_df, tibble(cls = train[[target_col]]))
  
  rule_tbl <- dat %>%
    count(across(starts_with("V")), cls, name = "n") %>%
    group_by(across(starts_with("V"))) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(out_mf = ifelse(cls == 0, 1L, 2L),
           weight = pmax(0.5, n / max(n)))
  
  rule_mat <- as.matrix(cbind(
    as.matrix(rule_tbl[, paste0("V", seq_along(feat_cols))]),
    rule_tbl$out_mf, rule_tbl$weight, 1L
  ))
  
  fis <- addrule(fis, rbind(rule_mat, add_default_rule(rule_mat, length(feat_cols))))
  
  list(fis = fis, n_rules_total = nrow(rule_mat) + 1)
}

n_mf_vec <- sapply(train[feat_cols], infer_nmf)
fls <- build_fls(train, feat_cols, target_col, n_mf_vec)

eval_model <- function(fis, newdata, feat_cols, y_true) {
  preds <- safe_eval(newdata[, feat_cols], fis)
  y_hat <- ifelse(preds > 0.5, 1, 0)
  list(conf = table(Predicted = y_hat, Actual = y_true),
       acc  = mean(y_hat == y_true))
}

res_fls <- eval_model(fls$fis, test, feat_cols, test[[target_col]])
cat("\n=== FLS (6-in) ===\n")
print(res_fls$conf)
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Rules:", fls$n_rules_total, "\n")

# -------------------------------
# 3) 2ISO builder (re-usable)
# -------------------------------
build_2iso <- function(df, in1, in2, target_col) {
  fis <- newfis(paste0("FIS_", in1, "_", in2), fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  
  feat_pair <- c(in1, in2)
  n_mf_vec <- sapply(df[feat_pair], infer_nmf)
  
  mf_defs <- vector("list", 2)
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
  
  fis <- addvar(fis, "output", "Risk", c(0,1))
  fis <- addmf(fis, "output", 1, "Benign",    "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Malignant", "trimf", c(0.5,1,1))
  
  idx_mat <- t(apply(as.matrix(df[, feat_pair]), 1, function(r) {
    sapply(seq_along(feat_pair), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  idx_df <- as.data.frame(idx_mat); colnames(idx_df) <- c("V1","V2")
  
  rule_tbl <- bind_cols(idx_df, tibble(cls = df[[target_col]])) %>%
    count(V1, V2, cls, name = "n") %>%
    group_by(V1, V2) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(out_mf = ifelse(cls == 0, 1L, 2L),
           weight = pmax(0.5, n / max(n)))
  
  rule_mat <- as.matrix(cbind(rule_tbl$V1, rule_tbl$V2,
                              rule_tbl$out_mf, rule_tbl$weight, 1L))
  fis <- addrule(fis, rbind(rule_mat, add_default_rule(rule_mat, 2)))
  
  list(fis = fis, n_rules_total = nrow(rule_mat) + 1)
}

# -------------------------------
# 4) HFS Serial (5 subsystems, 5 layers)
# -------------------------------
# S1: (radius, perimeter) -> size_risk
S1_ser <- build_2iso(train, "radius_mean", "perimeter_mean", target_col)
train$size_risk <- safe_eval(train[, c("radius_mean","perimeter_mean")], S1_ser$fis)
test$size_risk  <- safe_eval(test[,  c("radius_mean","perimeter_mean")], S1_ser$fis)

# S2: (size_risk, area) -> morphological_risk
S2_ser <- build_2iso(train, "size_risk", "area_mean", target_col)
train$morphological_risk <- safe_eval(train[, c("size_risk","area_mean")], S2_ser$fis)
test$morphological_risk  <- safe_eval(test[,  c("size_risk","area_mean")], S2_ser$fis)

# S3: (morphological_risk, texture) -> morphotexture_risk
S3_ser <- build_2iso(train, "morphological_risk", "texture_mean", target_col)
train$morphotexture_risk <- safe_eval(train[, c("morphological_risk","texture_mean")], S3_ser$fis)
test$morphotexture_risk  <- safe_eval(test[,  c("morphological_risk","texture_mean")], S3_ser$fis)

# S4: (morphotexture_risk, smoothness) -> surface_risk
S4_ser <- build_2iso(train, "morphotexture_risk", "smoothness_mean", target_col)
train$surface_risk <- safe_eval(train[, c("morphotexture_risk","smoothness_mean")], S4_ser$fis)
test$surface_risk  <- safe_eval(test[,  c("morphotexture_risk","smoothness_mean")], S4_ser$fis)

# S5: (surface_risk, compactness) -> Cancer Risk (final)
S5_ser <- build_2iso(train, "surface_risk", "compactness_mean", target_col)
pred_ser <- safe_eval(test[, c("surface_risk","compactness_mean")], S5_ser$fis)
yhat_ser <- ifelse(pred_ser > 0.5, 1, 0)
conf_ser <- table(Predicted = yhat_ser, Actual = test[[target_col]])
acc_ser  <- mean(yhat_ser == test[[target_col]])

# -------------------------------
# 5) HFS Parallel (5 subsystems, 4 layers)
# -------------------------------
# Layer 1
P1 <- build_2iso(train, "radius_mean", "perimeter_mean", target_col)   # size_risk
P2 <- build_2iso(train, "area_mean", "texture_mean", target_col)       # area_texture_risk
train$size_risk_p <- safe_eval(train[, c("radius_mean","perimeter_mean")], P1$fis)
test$size_risk_p  <- safe_eval(test[,  c("radius_mean","perimeter_mean")], P1$fis)
train$area_texture_risk <- safe_eval(train[, c("area_mean","texture_mean")], P2$fis)
test$area_texture_risk  <- safe_eval(test[,  c("area_mean","texture_mean")], P2$fis)

# Layer 2
P3 <- build_2iso(train, "smoothness_mean", "compactness_mean", target_col)  # surface_compactness_risk
train$surface_compactness_risk <- safe_eval(train[, c("smoothness_mean","compactness_mean")], P3$fis)
test$surface_compactness_risk  <- safe_eval(test[,  c("smoothness_mean","compactness_mean")], P3$fis)

# Layer 3
P4 <- build_2iso(train, "size_risk_p", "area_texture_risk", target_col)     # morphology_risk
train$morphology_risk <- safe_eval(train[, c("size_risk_p","area_texture_risk")], P4$fis)
test$morphology_risk  <- safe_eval(test[,  c("size_risk_p","area_texture_risk")], P4$fis)

# Layer 4
P5 <- build_2iso(train, "morphology_risk", "surface_compactness_risk", target_col)  # final
pred_par <- safe_eval(test[, c("morphology_risk","surface_compactness_risk")], P5$fis)
yhat_par <- ifelse(pred_par > 0.5, 1, 0)
conf_par <- table(Predicted = yhat_par, Actual = test[[target_col]])
acc_par  <- mean(yhat_par == test[[target_col]])

# -------------------------------
# 6) Summaries
# -------------------------------
cat("\n=== HFS Serial (5 subsystems, 5 layers) ===\n")
print(conf_ser)
cat("Accuracy:", sprintf("%.4f", acc_ser), "\n")
rules_ser <- S1_ser$n_rules_total + S2_ser$n_rules_total + S3_ser$n_rules_total +
  S4_ser$n_rules_total + S5_ser$n_rules_total
cat("Total Rules:", rules_ser, "\n")

cat("\n=== HFS Parallel (5 subsystems, 4 layers) ===\n")
print(conf_par)
cat("Accuracy:", sprintf("%.4f", acc_par), "\n")
rules_par <- P1$n_rules_total + P2$n_rules_total + P3$n_rules_total + P4$n_rules_total + P5$n_rules_total
cat("Total Rules:", rules_par, "\n")

cat("\n=== FLS (6-in) ===\n")
print(res_fls$conf)
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Total Rules:", fls$n_rules_total, "\n")

# LaTeX confusion matrices
cat("\nLaTeX — FLS:", to_latex_mat(res_fls$conf), "\n")
cat("LaTeX — HFS Serial:", to_latex_mat(conf_ser), "\n")
cat("LaTeX — HFS Parallel:", to_latex_mat(conf_par), "\n")

# -------------------------------
# 7) Performance table (LaTeX)
# -------------------------------
acc_fls <- res_fls$acc
acc_serial <- acc_ser
acc_parallel <- acc_par

latex_table <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Performance comparison between FLS, HFS Serial, and HFS Parallel for the Breast Cancer dataset (binary: Benign/Malignant). Six input variables were used: \\textit{radius\\_mean}, \\textit{texture\\_mean}, \\textit{perimeter\\_mean}, \\textit{area\\_mean}, \\textit{smoothness\\_mean}, \\textit{compactness\\_mean}.}\n",
  "\\label{tab:breast_results}\n",
  "\\resizebox{\\textwidth}{!}{\n",
  "\\begin{tabular}{l c c c c c c c c c}\n",
  "\\toprule\n",
  "\\textbf{Model} & \\textbf{No. of Inputs} & \\textbf{No. of Outputs} & \\textbf{No. of Subsystems} & \\textbf{No. of Layers} & \\textbf{Accuracy} & \\textbf{Confusion Matrix (Pred $\\times$ Actual)} & \\textbf{No. of Rules} & \\textbf{Avg. Rule Length} & \\textbf{Max Antecedents} \\\\\n",
  "\\midrule\n",
  sprintf("FLS (Data-driven) & 6 & 1 & 1 & 1 & %.3f & %s & %d & 6.0 & 6 \\\\\n",
          acc_fls, to_latex_mat(res_fls$conf), fls$n_rules_total),
  sprintf("HFS Serial & 6 & 1 & 5 & 5 & %.3f & %s & %d & 2.0 & 2 \\\\\n",
          acc_serial, to_latex_mat(conf_ser), rules_ser),
  sprintf("HFS Parallel & 6 & 1 & 5 & 4 & %.3f & %s & %d & 2.0 & 2 \\\\\n",
          acc_parallel, to_latex_mat(conf_par), rules_par),
  "\\bottomrule\n",
  "\\end{tabular}}\n",
  "\\end{table}\n"
)

cat("\n--- LaTeX Performance Table ---\n")
cat(latex_table)
cat("\n--------------------------------\n")

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
layers_fls <- list( list(fls$fis) )

# HFS Serial: 5 layers, 1 subsystem per layer
layers_ser <- list(
  list(S1_ser$fis),
  list(S2_ser$fis),
  list(S3_ser$fis),
  list(S4_ser$fis),
  list(S5_ser$fis)
)

# HFS Parallel: 4 layers (L1 has two subsystems)
layers_par <- list(
  list(P1$fis, P2$fis),
  list(P3$fis),
  list(P4$fis),
  list(P5$fis)
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
