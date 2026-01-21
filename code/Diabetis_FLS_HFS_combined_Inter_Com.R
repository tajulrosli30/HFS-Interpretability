#############################################
# BRFSS Diabetes (0:No, 1:Pre, 2:Diabetes)
# FLS (4-in), HFS-Parallel, HFS-Serial
# Stratified n=2000, Gaussian MFs via k-means,
# weighted rules + default rule
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
data_path  <- "G:/Other computers/My Laptop/UiTM - TAJUL/Postdoc - UoN/Publication/Paper HFS - Medical/Data/diabetes_012_health_indicators_BRFSS2023.csv"
target_col <- "diabetes_012"                       # 0,1,2
feat_cols  <- c("bmi","high_bp","high_chol","age_group")
set.seed(42)

df <- read_csv(data_path, show_col_types = FALSE) |>
  clean_names() |>
  select(all_of(c(feat_cols, target_col))) |>
  drop_na()

# -------------------------------
# Stratified sample of 2,000 rows
# -------------------------------
n_total <- 2000L
cls <- df[[target_col]]
tab <- as.integer(table(cls))
lev <- as.integer(names(table(cls)))
prop <- tab / sum(tab)

raw_targets <- prop * n_total
sizes <- floor(raw_targets)
rem <- n_total - sum(sizes)
if (rem > 0) {
  frac <- raw_targets - sizes
  bump_order <- order(frac, decreasing = TRUE)
  for (i in seq_len(rem)) sizes[bump_order[i]] <- sizes[bump_order[i]] + 1L
}
avail <- as.integer(table(cls))
sizes <- pmin(sizes, avail)

idx_keep <- integer(0)
for (j in seq_along(lev)) {
  idx_class <- which(df[[target_col]] == lev[j])
  take <- sizes[j]
  if (take > 0) idx_keep <- c(idx_keep, sample(idx_class, size = take, replace = FALSE))
}
df <- df[idx_keep, , drop = FALSE]
if (nrow(df) > n_total) df <- df[sample(seq_len(nrow(df)), n_total), , drop = FALSE]

# Stratified train/test 80/20
trainIndex <- createDataPartition(df[[target_col]], p = 0.8, list = FALSE)
train <- df[trainIndex, ]
test  <- df[-trainIndex, ]

# Ensure numeric inputs
train <- train %>% mutate(across(all_of(feat_cols), as.numeric))
test  <- test  %>% mutate(across(all_of(feat_cols), as.numeric))

# -------------------------------
# 1) Helpers: MFs & clustering
# -------------------------------
# Decide MF count per feature (binary-ish -> 2, else 3)
infer_nmf <- function(x) {
  ux <- sort(unique(x))
  if (length(ux) <= 3 || all(ux %in% c(0,1))) 2L else 3L
}
n_mf_vec <- sapply(train[feat_cols], infer_nmf)

# k-means centers (sorted), robust sigma
km_centers_sigma <- function(x, k) {
  x <- as.numeric(x)
  # guard for degenerate variance
  if (sd(x, na.rm = TRUE) == 0 || length(unique(x)) < k) {
    c_min <- min(x, na.rm = TRUE); c_max <- max(x, na.rm = TRUE)
    centers <- if (k == 2) c(c_min, c_max) else seq(c_min, c_max, length.out = k)
  } else {
    km <- kmeans(x, centers = k, nstart = 20)
    centers <- sort(as.numeric(km$centers))
  }
  # sigma: local neighborhood width, fallback to global range/(2*k)
  if (k == 1) {
    sigma <- max(sd(x, na.rm = TRUE), diff(range(x, na.rm = TRUE))/6, 1e-6)
    sigmas <- rep(sigma, 1)
  } else {
    gaps <- diff(centers)
    base_sigma <- if (all(gaps > 0)) {
      c(gaps[1], gaps[-1] + gaps[-length(gaps)], gaps[length(gaps)]) / 2
    } else {
      rep(diff(range(x, na.rm = TRUE))/(2*k), k)
    }
    # smooth & floor sigma
    sigmas <- pmax(base_sigma, diff(range(x, na.rm = TRUE))/(5*k))
  }
  list(centers = centers, sigmas = sigmas)
}

# Gaussian membership degree
gauss_eval <- function(x, sigma, c) exp(-0.5 * ((x - c)/sigma)^2)

# Build input MFs definition for one feature
build_feature_mfs <- function(x, k) {
  ini <- km_centers_sigma(x, k)
  centers <- ini$centers
  sigmas  <- ini$sigmas
  mfs <- vector("list", length = k)
  for (i in seq_len(k)) {
    mfs[[i]] <- list(type = "gaussmf", params = c(sigmas[i], centers[i]))
  }
  mfs
}

# Argmax MF index for a value
mf_argmax <- function(x, mf_list) {
  degs <- sapply(mf_list, function(m) gauss_eval(x, m$params[1], m$params[2]))
  which.max(degs)
}

# Map crisp output y in [0,2] to nearest class 0/1/2
crisp_to_class <- function(y) {
  classes <- c(0,1,2)
  classes[which.min(abs(y - classes))]
}

# Add default rule (low weight) if desired
add_default_rule <- function(rule_mat, n_in) {
  # antecedents = 1 for all inputs (i.e., first MF of each input)
  # but since we learn per-combination, default here means: a single low-weight rule per model
  # We'll implement as: no antecedent constraints, by repeating a common rule with low weight.
  # FuzzyR needs specific antecedents; we’ll use the first MF of each input by convention.
  r <- rep(1L, n_in)
  cbind(matrix(r, nrow = 1),  # antecedent indices
        output = 2L,          # class 1 (PreDiabetes) as neutral center
        weight = 0.1,         # low weight
        connector = 1L)
}

# -------------------------------
# 2) FLS (4 inputs -> Risk)
# -------------------------------
build_fls_4in <- function(train, feat_cols, target_col, n_mf_vec, add_default = TRUE) {
  # Build FIS shell with Gaussian MFs learned from data
  fis <- newfis("FLS_4in_3class", fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  
  # inputs
  mf_defs <- vector("list", length(feat_cols)); names(mf_defs) <- feat_cols
  for (i in seq_along(feat_cols)) {
    x <- train[[feat_cols[i]]]
    rng <- range(x, na.rm = TRUE)
    fis <- addvar(fis, "input", feat_cols[i], rng)
    mfs <- build_feature_mfs(x, n_mf_vec[i])
    mf_defs[[i]] <- mfs
    for (j in seq_along(mfs)) {
      par <- mfs[[j]]$params
      fis <- addmf(fis, "input", i, paste0("G", j), "gaussmf", par)
    }
  }
  
  # output (Mamdani): triangular around 0,1,2 on [0,2]
  fis <- addvar(fis, "output", "Risk", c(0,2))
  fis <- addmf(fis, "output", 1, "NoDiabetes",  "trimf", c(0, 0, 1))
  fis <- addmf(fis, "output", 1, "PreDiabetes", "trimf", c(0, 1, 2))
  fis <- addmf(fis, "output", 1, "Diabetes",    "trimf", c(1, 2, 2))
  
  # Data-driven weighted rules
  train_mat <- as.matrix(train[, feat_cols])
  idx_mat <- t(apply(train_mat, 1, function(r) {
    sapply(seq_along(feat_cols), function(i) mf_argmax(r[i], mf_defs[[i]]))
  }))
  colnames(idx_mat) <- paste0("i", seq_along(feat_cols))
  dat <- bind_cols(as_tibble(idx_mat), tibble(cls = train[[target_col]]))
  
  # Majority class per unique combo + weight = frequency
  rule_tbl <- dat |>
    count(across(starts_with("i")), cls, name = "n") |>
    group_by(across(starts_with("i"))) |>
    slice_max(order_by = n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(out_mf = case_when(cls == 0 ~ 1L, cls == 1 ~ 2L, TRUE ~ 3L),
           weight = pmax(0.5, n / max(n))) |>
    select(starts_with("i"), out_mf, weight) 
  
  rule_mat <- as.matrix(cbind(as.matrix(rule_tbl[, 1:length(feat_cols)]),
                              out_mf = rule_tbl$out_mf,
                              weight = rule_tbl$weight,
                              connector = 1L))
  storage.mode(rule_mat) <- "numeric"
  
  if (add_default) {
    def_rule <- add_default_rule(rule_mat, length(feat_cols))
    rule_mat <- rbind(rule_mat, def_rule)
  }
  
  fis <- addrule(fis, rule_mat)
  list(fis = fis,
       n_rules = nrow(rule_mat) - as.integer(add_default),
       n_rules_total = nrow(rule_mat),
       avg_rule_len = length(feat_cols),
       max_ante = length(feat_cols))
}

eval_model <- function(fis, newdata, feat_cols, y_true) {
  X <- as.matrix(newdata[, feat_cols])
  preds <- apply(X, 1, function(x) evalfis(x, fis))
  y_hat <- sapply(preds, crisp_to_class)
  conf <- table(Predicted = y_hat, Actual = y_true)
  acc <- mean(y_hat == y_true)
  list(conf = conf, acc = acc)
}

# -------------------------------
# 3) 2ISO builder & learner (used for HFS)
# -------------------------------
build_2iso_gauss <- function(name, x1_name, x2_name, x1_vec, x2_vec) {
  n1 <- infer_nmf(x1_vec); n2 <- infer_nmf(x2_vec)
  fis <- newfis(name, fisType = "mamdani",
                andMethod = "min", orMethod  = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  
  # x1
  rng1 <- range(x1_vec, na.rm = TRUE)
  fis <- addvar(fis, "input", x1_name, rng1)
  mfs1 <- build_feature_mfs(x1_vec, n1)
  for (j in seq_along(mfs1)) fis <- addmf(fis, "input", 1, paste0("G", j), "gaussmf", mfs1[[j]]$params)
  
  # x2
  rng2 <- range(x2_vec, na.rm = TRUE)
  fis <- addvar(fis, "input", x2_name, rng2)
  mfs2 <- build_feature_mfs(x2_vec, n2)
  for (j in seq_along(mfs2)) fis <- addmf(fis, "input", 2, paste0("G", j), "gaussmf", mfs2[[j]]$params)
  
  # output on [0,2]
  fis <- addvar(fis, "output", "Risk", c(0,2))
  fis <- addmf(fis, "output", 1, "NoDiabetes",  "trimf", c(0,0,1))
  fis <- addmf(fis, "output", 1, "PreDiabetes", "trimf", c(0,1,2))
  fis <- addmf(fis, "output", 1, "Diabetes",    "trimf", c(1,2,2))
  
  list(fis = fis, mfs = list(mfs1 = mfs1, mfs2 = mfs2), n1 = n1, n2 = n2)
}

learn_rules_2iso_weighted <- function(fis, mfs, df2, x1, x2, ycol, add_default = TRUE) {
  # fuzzify with argmax of Gaussian MFs
  idx_mat <- t(apply(as.matrix(df2[, c(x1, x2)]), 1, function(r) {
    c(mf_argmax(r[1], mfs$mfs1), mf_argmax(r[2], mfs$mfs2))
  }))
  colnames(idx_mat) <- c("i1","i2")
  dat <- bind_cols(as_tibble(idx_mat), tibble(cls = df2[[ycol]]))
  
  rule_tbl <- dat |>
    count(i1, i2, cls, name = "n") |>
    group_by(i1, i2) |>
    slice_max(order_by = n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(out_mf = case_when(cls == 0 ~ 1L, cls == 1 ~ 2L, TRUE ~ 3L),
           weight = pmax(0.5, n / max(n))) |>
    select(i1, i2, out_mf, weight)
  
  rule_mat <- as.matrix(cbind(rule_tbl$i1, rule_tbl$i2, rule_tbl$out_mf, rule_tbl$weight, 1L))
  storage.mode(rule_mat) <- "numeric"
  
  if (add_default) {
    def_rule <- c(1L, 1L, 2L, 0.1, 1L)
    rule_mat <- rbind(rule_mat, def_rule)
  }
  
  fis2 <- addrule(fis, rule_mat)
  list(fis = fis2, n_rules = nrow(rule_mat) - as.integer(add_default), n_rules_total = nrow(rule_mat))
}

# -------------------------------
# 4) Build & Evaluate MODELS
# -------------------------------

# ---- FLS 4-in
fls <- build_fls_4in(train, feat_cols, target_col, n_mf_vec, add_default = TRUE)
res_fls <- eval_model(fls$fis, test, feat_cols, test[[target_col]])

cat("\n=== FLS (4-in) ===\n")
print(res_fls$conf)
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Rules:", fls$n_rules, "(+ defaults =", fls$n_rules_total - fls$n_rules, ")",
    "| Total:", fls$n_rules_total, "\n")
cat("Avg. Rule Length:", fls$avg_rule_len, "| Max Antecedents:", fls$max_ante, "\n")

# ---- HFS Parallel
# Layer1: S1(bmi,high_bp)  S2(high_chol,age_group)
# Layer2: S3(y1,y2)
S1 <- build_2iso_gauss("S1_bmi_bp", "bmi","high_bp", train$bmi, train$high_bp)
S1_learn <- learn_rules_2iso_weighted(S1$fis, S1$mfs, train, "bmi","high_bp", target_col)

S2 <- build_2iso_gauss("S2_chol_age", "high_chol","age_group", train$high_chol, train$age_group)
S2_learn <- learn_rules_2iso_weighted(S2$fis, S2$mfs, train, "high_chol","age_group", target_col)

# get y1,y2 on TRAIN
y1_tr <- apply(as.matrix(train[, c("bmi","high_bp")]), 1, function(x) evalfis(x, S1_learn$fis))
y2_tr <- apply(as.matrix(train[, c("high_chol","age_group")]), 1, function(x) evalfis(x, S2_learn$fis))

# S3 with y1,y2 (range [0,2]); make simple 3-Gaussian per input anchored at 0,1,2
mk_y_mfs <- function(fis_name) {
  fis <- newfis(fis_name, fisType = "mamdani",
                andMethod = "min", orMethod = "max",
                impMethod = "min", aggMethod = "max",
                defuzzMethod = "centroid")
  fis <- addvar(fis, "input", "y1", c(0,2))
  for (mu in c(0,1,2)) fis <- addmf(fis, "input", 1, paste0("G", mu), "gaussmf", c(0.35, mu))
  fis <- addvar(fis, "input", "y2", c(0,2))
  for (mu in c(0,1,2)) fis <- addmf(fis, "input", 2, paste0("G", mu), "gaussmf", c(0.35, mu))
  fis <- addvar(fis, "output", "Risk", c(0,2))
  fis <- addmf(fis, "output", 1, "NoDiabetes",  "trimf", c(0,0,1))
  fis <- addmf(fis, "output", 1, "PreDiabetes", "trimf", c(0,1,2))
  fis <- addmf(fis, "output", 1, "Diabetes",    "trimf", c(1,2,2))
  fis
}

S3 <- mk_y_mfs("S3_y1y2")
# Learn S3 rules (y1,y2)->class using coarse MF argmax (0/1/2)
argmax_012 <- function(v) which.min(abs(c(0,1,2) - v)) # returns 1..3
idx_y <- cbind(i1 = sapply(y1_tr, argmax_012),
               i2 = sapply(y2_tr, argmax_012))
dat_y <- bind_cols(as_tibble(idx_y), tibble(cls = train[[target_col]]))
tbl_y <- dat_y |>
  count(i1, i2, cls, name = "n") |>
  group_by(i1, i2) |>
  slice_max(order_by = n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(out_mf = case_when(cls == 0 ~ 1L, cls == 1 ~ 2L, TRUE ~ 3L),
         weight = pmax(0.5, n / max(n))) |>
  select(i1, i2, out_mf, weight)
rule_S3 <- as.matrix(cbind(tbl_y$i1, tbl_y$i2, tbl_y$out_mf, tbl_y$weight, 1L))
rule_S3 <- rbind(rule_S3, c(2L,2L,2L,0.1,1L)) # default to center
S3 <- addrule(S3, rule_S3)

# Evaluate HFS-Parallel on TEST
y1_te <- apply(as.matrix(test[, c("bmi","high_bp")]), 1, function(x) evalfis(x, S1_learn$fis))
y2_te <- apply(as.matrix(test[, c("high_chol","age_group")]), 1, function(x) evalfis(x, S2_learn$fis))
risk_par <- mapply(function(a,b) evalfis(c(a,b), S3), y1_te, y2_te)
yhat_par <- sapply(risk_par, crisp_to_class)
conf_par <- table(Predicted = yhat_par, Actual = test[[target_col]])
acc_par  <- mean(yhat_par == test[[target_col]])

cat("\n=== HFS Parallel ===\n")
print(conf_par)
cat("Accuracy:", sprintf("%.4f", acc_par), "\n")
rules_par <- c(S1 = S1_learn$n_rules_total, S2 = S2_learn$n_rules_total, S3 = nrow(rule_S3))
cat("Rules per subsystem:", paste(names(rules_par), rules_par, collapse = " | "),
    "| Total:", sum(rules_par), "\n")
cat("Avg. Rule Length: 2 | Max Antecedents: 2 | Subsystems: 3 | Layers: 2\n")

# ---- HFS Serial
# y1 = SS1(bmi, high_bp)
SS1 <- build_2iso_gauss("SS1_bmi_bp", "bmi","high_bp", train$bmi, train$high_bp)
SS1_learn <- learn_rules_2iso_weighted(SS1$fis, SS1$mfs, train, "bmi","high_bp", target_col)

# y2 = SS2(y1, high_chol)
y1_tr_s <- apply(as.matrix(train[, c("bmi","high_bp")]), 1, function(x) evalfis(x, SS1_learn$fis))
SS2 <- mk_y_mfs("SS2_y1chol") # reuse structure but rename inputs to y1 & high_chol
# rebuild SS2 with (y1 in [0,2], high_chol from data)
SS2 <- newfis("SS2_y1chol", fisType = "mamdani",
              andMethod = "min", orMethod  = "max",
              impMethod = "min", aggMethod = "max",
              defuzzMethod = "centroid")
# y1
SS2 <- addvar(SS2, "input", "y1", c(0,2))
for (mu in c(0,1,2)) SS2 <- addmf(SS2, "input", 1, paste0("G", mu), "gaussmf", c(0.35, mu))
# high_chol
hc_rng <- range(train$high_chol, na.rm = TRUE)
SS2 <- addvar(SS2, "input", "high_chol", hc_rng)
# use inferred MF count for high_chol
n_hc <- infer_nmf(train$high_chol)
mfs_hc <- build_feature_mfs(train$high_chol, n_hc)
for (j in seq_along(mfs_hc)) SS2 <- addmf(SS2, "input", 2, paste0("G", j), "gaussmf", mfs_hc[[j]]$params)
# output
SS2 <- addvar(SS2, "output", "Risk", c(0,2))
SS2 <- addmf(SS2, "output", 1, "NoDiabetes",  "trimf", c(0,0,1))
SS2 <- addmf(SS2, "output", 1, "PreDiabetes", "trimf", c(0,1,2))
SS2 <- addmf(SS2, "output", 1, "Diabetes",    "trimf", c(1,2,2))
# learn SS2 rules
train_SS2 <- train %>% mutate(y1 = y1_tr_s)

# fuzzify (y1 via 0/1/2 argmax; high_chol via mfs_hc)
idx_SS2 <- t(apply(as.matrix(train_SS2[, c("y1","high_chol")]), 1, function(r) {
  c(argmax_012(r[1]), mf_argmax(r[2], mfs_hc))
}))

# Give proper names so count() works
colnames(idx_SS2) <- c("V1", "V2")

dat_SS2 <- bind_cols(as_tibble(idx_SS2), tibble(cls = train[[target_col]]))

tbl_SS2 <- dat_SS2 |>
  count(V1, V2, cls, name = "n") |>
  group_by(V1, V2) |>
  slice_max(order_by = n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(out_mf = case_when(cls == 0 ~ 1L, cls == 1 ~ 2L, TRUE ~ 3L),
         weight = pmax(0.5, n / max(n)))
rule_SS2 <- as.matrix(cbind(tbl_SS2$V1, tbl_SS2$V2, tbl_SS2$out_mf, tbl_SS2$weight, 1L))
rule_SS2 <- rbind(rule_SS2, c(2L,1L,2L,0.1,1L))
SS2 <- addrule(SS2, rule_SS2)

# final = SS3(y2, age_group)
y2_tr_s <- mapply(function(a,b) evalfis(c(a,b), SS2), y1_tr_s, train$high_chol)
SS3 <- newfis("SS3_y2age", fisType = "mamdani",
              andMethod = "min", orMethod  = "max",
              impMethod = "min", aggMethod = "max",
              defuzzMethod = "centroid")
# y2
SS3 <- addvar(SS3, "input", "y2", c(0,2))
for (mu in c(0,1,2)) SS3 <- addmf(SS3, "input", 1, paste0("G", mu), "gaussmf", c(0.35, mu))
# age_group
ag_rng <- range(train$age_group, na.rm = TRUE)
SS3 <- addvar(SS3, "input", "age_group", ag_rng)
n_ag <- infer_nmf(train$age_group)
mfs_ag <- build_feature_mfs(train$age_group, n_ag)
for (j in seq_along(mfs_ag)) SS3 <- addmf(SS3, "input", 2, paste0("G", j), "gaussmf", mfs_ag[[j]]$params)
# output
SS3 <- addvar(SS3, "output", "Risk", c(0,2))
SS3 <- addmf(SS3, "output", 1, "NoDiabetes",  "trimf", c(0,0,1))
SS3 <- addmf(SS3, "output", 1, "PreDiabetes", "trimf", c(0,1,2))
SS3 <- addmf(SS3, "output", 1, "Diabetes",    "trimf", c(1,2,2))
# fuzzify (y2 via argmax; age_group via mfs_age)
train_SS3 <- train %>% mutate(y2 = y2_tr_s)  # <-- you also need to define train_SS3 here
idx_SS3 <- t(apply(as.matrix(train_SS3[, c("y2","age_group")]), 1, function(r) {
  c(argmax_012(r[1]), mf_argmax(r[2], mfs_ag))
}))
# Ensure proper column names
colnames(idx_SS3) <- c("V1", "V2")

dat_SS3 <- bind_cols(as_tibble(idx_SS3), tibble(cls = train[[target_col]]))

tbl_SS3 <- dat_SS3 |>
  count(V1, V2, cls, name = "n") |>
  group_by(V1, V2) |>
  slice_max(order_by = n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(out_mf = case_when(cls == 0 ~ 1L, cls == 1 ~ 2L, TRUE ~ 3L),
         weight = pmax(0.5, n / max(n)))
rule_SS3 <- as.matrix(cbind(tbl_SS3$V1, tbl_SS3$V2, tbl_SS3$out_mf, tbl_SS3$weight, 1L))
rule_SS3 <- rbind(rule_SS3, c(2L,1L,2L,0.1,1L))
SS3 <- addrule(SS3, rule_SS3)

# Evaluate HFS-Serial on TEST
y1_te_s <- apply(as.matrix(test[, c("bmi","high_bp")]), 1, function(x) evalfis(x, SS1_learn$fis))
y2_te_s <- mapply(function(a,b) evalfis(c(a,b), SS2), y1_te_s, test$high_chol)
risk_ser <- mapply(function(a,b) evalfis(c(a,b), SS3), y2_te_s, test$age_group)
yhat_ser <- sapply(risk_ser, crisp_to_class)
conf_ser <- table(Predicted = yhat_ser, Actual = test[[target_col]])
acc_ser  <- mean(yhat_ser == test[[target_col]])

cat("\n=== HFS Serial ===\n")
print(conf_ser)
cat("Accuracy:", sprintf("%.4f", acc_ser), "\n")
rules_ser <- c(SS1 = nrow(rule_SS2) + 0*0 + 0,  # placeholder, compute correct
               SS2 = nrow(rule_SS2),
               SS3 = nrow(rule_SS3))
# Recompute SS1 rules properly:
rules_ser["SS1"] <- length(SS1_learn$n_rules_total)
cat("Rules per subsystem: SS1", SS1_learn$n_rules_total, 
    "| SS2", nrow(rule_SS2), "| SS3", nrow(rule_SS3),
    "| Total:", (SS1_learn$n_rules_total + nrow(rule_SS2) + nrow(rule_SS3)), "\n")
cat("Avg. Rule Length: 2 | Max Antecedents: 2 | Subsystems: 3 | Layers: 3\n")

# -------------------------------
# 5) Quick summary
# -------------------------------


cat("\n--- Summary ---\n")

# --- FLS (4-in)
cat("\nFLS (4-in)\n")
cat("Accuracy:", sprintf("%.4f", res_fls$acc), "\n")
cat("Rules per subsystem: Full system =", fls$n_rules_total, "\n")
cat("Total Rules:", fls$n_rules_total, "\n")
cat("Avg. Rule Length:", fls$avg_rule_len, 
    "| Max Antecedents:", fls$max_ante, 
    "| Subsystems: 1 | Layers: 1\n")

# --- HFS Parallel
cat("\nHFS Parallel\n")
cat("Accuracy:", sprintf("%.4f", acc_par), "\n")
cat("Rules per subsystem:", paste(names(rules_par), rules_par, collapse = " | "),
    "| Total:", sum(rules_par), "\n")
cat("Avg. Rule Length: 2 | Max Antecedents: 2 | Subsystems: 3 | Layers: 2\n")

# --- HFS Serial
cat("\nHFS Serial\n")
cat("Accuracy:", sprintf("%.4f", acc_ser), "\n")
rules_ser <- c(SS1 = SS1_learn$n_rules_total, 
               SS2 = nrow(rule_SS2), 
               SS3 = nrow(rule_SS3))
cat("Rules per subsystem:", paste(names(rules_ser), rules_ser, collapse = " | "),
    "| Total:", sum(rules_ser), "\n")
cat("Avg. Rule Length: 2 | Max Antecedents: 2 | Subsystems: 3 | Layers: 3\n")

# ---- Package HFS-Parallel results like 'res_hfs_parallel'
res_hfs_parallel <- list(
  acc = acc_par,
  conf = conf_par,
  n_rules_total = (S1_learn$n_rules_total + S2_learn$n_rules_total + nrow(rule_S3)),
  avg_rule_len = 2L,
  max_ante = 2L,
  n_subsystems = 3L,
  n_layers = 2L
)

# ---- Package HFS-Serial results like 'res_hfs_serial'
# (Fix: don't use length(); use the actual total)
rules_ser_vec <- c(
  SS1 = SS1_learn$n_rules_total,
  SS2 = nrow(rule_SS2),
  SS3 = nrow(rule_SS3)
)
res_hfs_serial <- list(
  acc = acc_ser,
  conf = conf_ser,
  n_rules_total = sum(rules_ser_vec),
  avg_rule_len = 2L,
  max_ante = 2L,
  n_subsystems = 3L,
  n_layers = 3L
)

# ---- Helpers to produce LaTeX matrix text
to_latex_mat <- function(tab) {
  m <- as.matrix(tab)
  apply(m, 1, function(r) paste(r, collapse = " & ")) |>
    (\(rows) paste0("\\begin{bmatrix}", paste(rows, collapse=" \\\\ "), " \\end{bmatrix}"))()
}

cat("\nHFS Parallel — Accuracy:", sprintf("%.4f", res_hfs_parallel$acc), "\n")
cat("Confusion (LaTeX): ", to_latex_mat(res_hfs_parallel$conf), "\n")
cat("Total Rules:", res_hfs_parallel$n_rules_total, "\n")

cat("\nHFS Serial — Accuracy:", sprintf("%.4f", res_hfs_serial$acc), "\n")
cat("Confusion (LaTeX): ", to_latex_mat(res_hfs_serial$conf), "\n")
cat("Total Rules:", res_hfs_serial$n_rules_total, "\n")

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

# HFS Serial: 3 layers, 1 subsystem per layer
layers_ser <- list(
  list(S1_ser$fis),
  list(S2_ser$fis),
  list(S3_ser$fis)
)

# HFS Parallel: first layer has 2 subsystems, second layer has 1
layers_par <- list(
  list(P1$fis, P2$fis),  # Layer 1
  list(P3$fis)           # Layer 2
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
