#############################################
# Liver Disease — FLS, HFS Parallel, HFS Serial
# Binary target: Dataset (1 = Liver disease, 2 = No liver disease)
# Inputs (7): Age, Total_Bilirubin, Direct_Bilirubin,
#             Alkaline_Phosphotase, Alamine_Aminotransferase,
#             Aspartate_Aminotransferase, Total_Protiens
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
data_path <- "G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/Data/indian_liver_patient.csv"
target_col <- "dataset"
feat_cols  <- c("age", "total_bilirubin", "direct_bilirubin",
                "alkaline_phosphotase", "alamine_aminotransferase",
                "aspartate_aminotransferase", "total_protiens")
set.seed(42)

# Read and clean names
df <- read_csv(data_path, show_col_types = FALSE) |> clean_names()

# Check actual names from CSV and rename if necessary
actual_names <- names(df)
name_map <- c(
  "age"                          = "age",
  "total_bilirubin"              = "total_bilirubin",
  "direct_bilirubin"             = "direct_bilirubin",
  "alkaline_phosphotase"         = "alkaline_phosphotase",
  "alamine_aminotransferase"     = "alamine_aminotransferase",
  "aspartate_aminotransferase"   = "aspartate_aminotransferase",
  "total_protiens"               = "total_protiens",
  "dataset"                      = "dataset"
)
df <- df |> rename(any_of(name_map)) |> select(all_of(c(feat_cols, target_col))) |> drop_na()

# Encode target: 1 = liver disease (positive), 2 = no disease (negative)
df <- df |> mutate(!!target_col := ifelse(.data[[target_col]] == 1, 1L, 0L))

# Balance dataset to help HFS performance
df_bal <- upSample(x = df[, feat_cols], y = as.factor(df[[target_col]]))
names(df_bal)[ncol(df_bal)] <- target_col
df_bal[[target_col]] <- as.integer(as.character(df_bal[[target_col]]))

# Stratified train/test split
trainIndex <- createDataPartition(df_bal[[target_col]], p = 0.8, list = FALSE)
train <- df_bal[trainIndex, ]
test  <- df_bal[-trainIndex, ]

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
  x <- x[!is.na(x)]
  if (length(unique(x)) < k) {
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
# 2) FLS (7-in) — data-driven rules
# -------------------------------
build_fls <- function(train, feat_cols, target_col, n_mf_vec) {
  fis <- newfis("FLS_7in", fisType = "mamdani",
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
  
  fis <- addvar(fis, "output", "DiseaseRisk", c(0,1))
  fis <- addmf(fis, "output", 1, "NoDisease", "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Disease",  "trimf", c(0.5,1,1))
  
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
  fis <- addmf(fis, "output", 1, "NoDisease", "trimf", c(0,0,0.5))
  fis <- addmf(fis, "output", 1, "Disease",  "trimf", c(0.5,1,1))
  
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
# 4) HFS Serial (6 subsystems, 6 layers)
# -------------------------------
S1_ser <- build_2iso(train, "age", "total_bilirubin", target_col)
train$r1 <- safe_eval(train[, c("age","total_bilirubin")], S1_ser$fis)
test$r1  <- safe_eval(test[,  c("age","total_bilirubin")], S1_ser$fis)

S2_ser <- build_2iso(train, "r1", "direct_bilirubin", target_col)
train$r2 <- safe_eval(train[, c("r1","direct_bilirubin")], S2_ser$fis)
test$r2  <- safe_eval(test[,  c("r1","direct_bilirubin")], S2_ser$fis)

S3_ser <- build_2iso(train, "r2", "alkaline_phosphotase", target_col)
train$r3 <- safe_eval(train[, c("r2","alkaline_phosphotase")], S3_ser$fis)
test$r3  <- safe_eval(test[,  c("r2","alkaline_phosphotase")], S3_ser$fis)

S4_ser <- build_2iso(train, "r3", "alamine_aminotransferase", target_col)
train$r4 <- safe_eval(train[, c("r3","alamine_aminotransferase")], S4_ser$fis)
test$r4  <- safe_eval(test[,  c("r3","alamine_aminotransferase")], S4_ser$fis)

S5_ser <- build_2iso(train, "r4", "aspartate_aminotransferase", target_col)
train$r5 <- safe_eval(train[, c("r4","aspartate_aminotransferase")], S5_ser$fis)
test$r5  <- safe_eval(test[,  c("r4","aspartate_aminotransferase")], S5_ser$fis)

S6_ser <- build_2iso(train, "r5", "total_protiens", target_col)
pred_ser <- safe_eval(test[, c("r5","total_protiens")], S6_ser$fis)
yhat_ser <- ifelse(pred_ser > 0.5, 1, 0)
conf_ser <- table(Predicted = yhat_ser, Actual = test[[target_col]])
acc_ser  <- mean(yhat_ser == test[[target_col]])
rules_ser <- S1_ser$n_rules_total + S2_ser$n_rules_total + S3_ser$n_rules_total +
  S4_ser$n_rules_total + S5_ser$n_rules_total + S6_ser$n_rules_total

# -------------------------------
# 5) HFS Parallel (6 subsystems, 5 layers)
# -------------------------------
P1 <- build_2iso(train, "age", "total_bilirubin", target_col)
train$p1 <- safe_eval(train[, c("age","total_bilirubin")], P1$fis)
test$p1  <- safe_eval(test[,  c("age","total_bilirubin")], P1$fis)

P2 <- build_2iso(train, "direct_bilirubin", "alkaline_phosphotase", target_col)
train$p2 <- safe_eval(train[, c("direct_bilirubin","alkaline_phosphotase")], P2$fis)
test$p2  <- safe_eval(test[,  c("direct_bilirubin","alkaline_phosphotase")], P2$fis)

P3 <- build_2iso(train, "alamine_aminotransferase", "aspartate_aminotransferase", target_col)
train$p3 <- safe_eval(train[, c("alamine_aminotransferase","aspartate_aminotransferase")], P3$fis)
test$p3  <- safe_eval(test[,  c("alamine_aminotransferase","aspartate_aminotransferase")], P3$fis)

P4 <- build_2iso(train, "p1", "p2", target_col)
train$p4 <- safe_eval(train[, c("p1","p2")], P4$fis)
test$p4  <- safe_eval(test[,  c("p1","p2")], P4$fis)

P5 <- build_2iso(train, "p3", "total_protiens", target_col)
train$p5 <- safe_eval(train[, c("p3","total_protiens")], P5$fis)
test$p5  <- safe_eval(test[,  c("p3","total_protiens")], P5$fis)

P6 <- build_2iso(train, "p4", "p5", target_col)
pred_par <- safe_eval(test[, c("p4","p5")], P6$fis)
yhat_par <- ifelse(pred_par > 0.5, 1, 0)
conf_par <- table(Predicted = yhat_par, Actual = test[[target_col]])
acc_par  <- mean(yhat_par == test[[target_col]])
rules_par <- P1$n_rules_total + P2$n_rules_total + P3$n_rules_total + P4$n_rules_total + P5$n_rules_total + P6$n_rules_total

# -------------------------------
# 6) Results
# -------------------------------
cat("\n=== HFS Serial (6 subsystems, 6 layers) ===\n")
print(conf_ser)
cat("Accuracy:", sprintf("%.4f", acc_ser), "\n")
cat("Total Rules:", rules_ser, "\n")

cat("\n=== HFS Parallel (6 subsystems, 5 layers) ===\n")
print(conf_par)
cat("Accuracy:", sprintf("%.4f", acc_par), "\n")
cat("Total Rules:", rules_par, "\n")

cat("\n=== FLS (7-in) ===\n")
print(res_fls$conf)
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Total Rules:", fls$n_rules_total, "\n")

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

# HFS Serial: 6 layers, 1 subsystem each (S1_ser ... S6_ser)
layers_ser <- list(
  list(S1_ser$fis),
  list(S2_ser$fis),
  list(S3_ser$fis),
  list(S4_ser$fis),
  list(S5_ser$fis),
  list(S6_ser$fis)
)

## ---- HFS Parallel options ----

# Option B (paper-aligned, 5 layers) — if you want "6 subsystems, 5 layers"
layers_par <- list(
  list(P1$fis, P2$fis),          # L1
  list(P3$fis),                  # L2
  list(P4$fis),                  # L3
  list(P5$fis),                  # L4
  list(P6$fis)                   # L5
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

