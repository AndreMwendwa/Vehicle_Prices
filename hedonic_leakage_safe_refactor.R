
# ==============================
# Hedonic Regression (Leakage-Safe)
# ==============================
# This script refactors the original code to eliminate data leakage in:
# - Clustering (make-model clusters and brand clusters)
# - Cluster-level and make/brand make-model average price features (lagged)
# It also keeps PCA training strictly on the training split.
#
# Key design:
# 1) Time split FIRST (rsample::initial_time_split).
# 2) All transforms (k-means, scaling, averages, PCA) fit on TRAIN ONLY.
# 3) Test set is only transformed/assigned using train-fitted artifacts.
#
# NOTE: Paths/column names follow your original script.
#       You may need to adjust file paths if your working directory differs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(FactoMineR)
  library(factoextra)
  library(multiColl)
  library(rsample)
  library(Metrics)
  library(openxlsx)
  library(lubridate)
})

set.seed(123)

# ------------------------------
# 0) Load & Preprocess Base Data
# ------------------------------
df <- read_csv("delistings_secondary_data.csv")

df <- df |>
  separate_wider_position(yr_month, c(yr = 4, mnth = 2)) |>
  mutate(
    yr_month_dtype = make_date(year = yr, month = mnth),
    yr = as.integer(yr),
    mnth = as.integer(mnth)
  ) |>
  arrange(yr_month_dtype)

# Observational time markers (keep if needed)
df <- df |>
  mutate(
    yr_obsvn_before_june_2022 = pmax(0, 2022 + 6/12 - yr - mnth/12),
    yr_obsvn_after_june_2022  = pmax(0, yr + mnth/12 - 2022 - 6/12),
    yr_decimal = yr + mnth/12
  )

# Compute age (requires a 'year' column as Date in your base file)
df <- df |>
  mutate(year = make_date(year = year),
         age  = as.numeric(time_length(yr_month_dtype - year, "years")))

# Factors and reference levels
df <- df |>
  mutate(
    mnth     = as.factor(mnth),
    province = as.factor(province),
    make     = as.factor(make),
    Body_type= as.factor(Body_type)
  )

df$make     <- fct_relevel(df$make, 'Toyota')
df$province <- fct_relevel(df$province, 'Ontario')
df$Body_type<- fct_relevel(df$Body_type, 'Sedan')

# ------------------------------
# 1) CPI Merge & Deflation
# ------------------------------
cpi <- read_csv('CPI_data/18100006.csv')

cpi_canada_avg <- cpi |>
  filter(GEO == 'Canada') |>
  mutate(REF_DATE_MOD = as.Date(paste0(REF_DATE, "-01"))) |>
  group_by(REF_DATE_MOD) |>
  summarise(avg_value = mean(VALUE, na.rm = TRUE), .groups = "drop")

df <- df |>
  inner_join(cpi_canada_avg, join_by(yr_month_dtype == REF_DATE_MOD))

# Use a fixed reference (example: the last available CPI in your original code)
cpi_current <- tail(cpi_canada_avg$avg_value, 5)[[1]]

df <- df |>
  mutate(prices_deflated = (price_10k_bucket + 0.5) * 10000 * cpi_current / avg_value)

# Winsorize some columns at 99th percentile
list_cols_created <- c("Dimensions_Length", "Dimensions_Width", "Dimensions_Wheelbase",
                       "Cargo_Volume", "Fuel_Economy", "Engine_Torque_Cleaned", "Acceleration")
for (col in list_cols_created) {
  q <- quantile(df[[col]], 0.99, na.rm = TRUE)
  df[[col]] <- ifelse(df[[col]] > q, q, df[[col]])
}

# ------------------------------
# 2) Train/Test Split FIRST
# ------------------------------
# Keep brand filtering (n > 100) **per condition** but do it pre-split per your logic
df_new  <- df |> filter(condition == "New")
df_used <- df |> filter(condition == "Used")

count_filter_new  <- df_new  |> count(make) |> filter(n > 100) |> pull(make)
count_filter_used <- df_used |> count(make) |> filter(n > 100) |> pull(make)

df_new  <- df_new  |> filter(make %in% count_filter_new)
df_used <- df_used |> filter(make %in% count_filter_used)

split_new  <- initial_time_split(df_new,  prop = 0.9)
split_used <- initial_time_split(df_used, prop = 0.9)

df_new_train  <- training(split_new)
df_new_test   <- testing(split_new)
df_used_train <- training(split_used)
df_used_test  <- testing(split_used)

# Test set basic cleaning (as in original)
df_new_test  <- df_new_test  |> filter(!is.na(Body_type), !is.na(province))
df_used_test <- df_used_test |> filter(!is.na(Body_type), !is.na(province))

# Make-model keys AFTER split
df_new_train  <- df_new_train  |> mutate(make_model  = paste0(tolower(make), tolower(model)))
df_new_test   <- df_new_test   |> mutate(make_model  = paste0(tolower(make), tolower(model)))
df_used_train <- df_used_train |> mutate(make_model = paste0(tolower(make), tolower(model)))
df_used_test  <- df_used_test  |> mutate(make_model = paste0(tolower(make), tolower(model)))

# ------------------------------
# 3) Helpers: Aggregation, KMeans on TRAIN, Assign to TEST
# ------------------------------
agg_specs_by <- function(df, by_col) {
  df |>
    dplyr::select(
      Engine_Displacement_Cleaned, Engine_Torque_Cleaned, Engine_Capacity, Top_Speed,
      Acceleration, Dimensions_Length, Dimensions_Width, Dimensions_Wheelbase,
      Cargo_Volume, Unladen_Weight, Fuel_Economy, !!rlang::sym(by_col)
    ) |>
    group_by(!!rlang::sym(by_col)) |>
    summarise(
      n = n(),
      across(
        c(Engine_Displacement_Cleaned, Engine_Torque_Cleaned, Engine_Capacity, Top_Speed,
          Acceleration, Dimensions_Length, Dimensions_Width, Dimensions_Wheelbase,
          Cargo_Volume, Unladen_Weight, Fuel_Economy),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )
}

fit_clusters_train <- function(df_train, by_col = "make_model", k = 7) {
  agg <- agg_specs_by(df_train, by_col)
  rn  <- agg[[by_col]]
  
  mat <- as.matrix(agg |> dplyr::select(-all_of(c(by_col, "n"))))
  means <- colMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 2, sd)
  sds[sds == 0] <- 1
  mat_sc <- scale(mat, center = means, scale = sds)
  
  set.seed(123)
  km <- kmeans(mat_sc, centers = k, nstart = 25)
  
  tibble_map <- tibble(
    !!by_col := rn,
    cluster = factor(km$cluster)
  )
  
  list(
    by_col = by_col,
    centers = km$centers,
    means   = means,
    sds     = sds,
    mapping = tibble_map
  )
}

assign_clusters <- function(df, model) {
  agg <- agg_specs_by(df, model$by_col)
  rn  <- agg[[model$by_col]]
  mat <- as.matrix(agg |> dplyr::select(-all_of(c(model$by_col, "n"))))

  mat_sc <- scale(mat, center = model$means, scale = model$sds)

  dist_to_centers <- function(x, centers) apply(centers, 1, function(cn) sqrt(sum((x - cn)^2)))
  cl_idx <- apply(mat_sc, 1, function(row) which.min(dist_to_centers(row, model$centers)))

  tibble(!!model$by_col := rn, cluster = factor(cl_idx))
}

# ------------------------------
# 4) Leakage-safe Average Price Features (TRAIN-only), Lagged
# ------------------------------
# Helper to compute (t-12) lag on TRAIN ONLY, then left_join to train & test
lagged_avg_by_group_train <- function(df_train, group_cols, price_col = "prices_deflated",
                                      date_col = "yr_month_dtype",
                                      lag_months = 12,
                                      output_name_current = "avg_price_current",
                                      output_name_lag = "avg_price_prev") {

  df_avg <- df_train |>
    group_by(across(all_of(c(date_col, group_cols)))) |>
    summarise("{output_name_current}" := mean(.data[[price_col]], na.rm = TRUE),
              .groups = "drop_last") |>
    arrange(across(all_of(c(group_cols, date_col)))) |>
    group_by(across(all_of(group_cols))) |>
    mutate("{output_name_lag}" := dplyr::lag(.data[[output_name_current]], 12)) |>
    ungroup()

  df_avg |>
    dplyr::select(all_of(c(date_col, group_cols, output_name_lag)))
}

# ------------------------------
# 5) USED: Make-Model Clusters (k=7), Train-only; Lagged Cluster Prices
# ------------------------------
used_clust <- fit_clusters_train(df_used_train, by_col = "make_model", k = 7)

df_used_train <- df_used_train |> left_join(used_clust$mapping, by = "make_model")
# assign on test by nearest centroid (no refit)
df_used_test  <- df_used_test  |>
  left_join(assign_clusters(df_used_test, used_clust), by = "make_model")

# Train-only cluster lag features
used_cluster_lag <- lagged_avg_by_group_train(
  df_train   = df_used_train,
  group_cols = "cluster",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_cluster",
  output_name_lag     = "avg_price_previous_period_cluster"
)

df_used_train <- df_used_train |>
  left_join(used_cluster_lag, by = c("yr_month_dtype", "cluster"))
df_used_test <- df_used_test |>
  left_join(used_cluster_lag, by = c("yr_month_dtype", "cluster"))

# ------------------------------
# 6) NEW: Make-Model Clusters (k=6), Train-only; Lagged Cluster Prices
# ------------------------------
new_clust <- fit_clusters_train(df_new_train, by_col = "make_model", k = 6)

df_new_train <- df_new_train |> left_join(new_clust$mapping, by = "make_model")
df_new_test  <- df_new_test  |>
  left_join(assign_clusters(df_new_test, new_clust), by = "make_model")

new_cluster_lag <- lagged_avg_by_group_train(
  df_train   = df_new_train,
  group_cols = "cluster",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_cluster",
  output_name_lag     = "avg_price_previous_period_cluster"
)

df_new_train <- df_new_train |>
  left_join(new_cluster_lag, by = c("yr_month_dtype", "cluster"))
df_new_test <- df_new_test |>
  left_join(new_cluster_lag, by = c("yr_month_dtype", "cluster"))

# ------------------------------
# 7) Brand-only Clusters (cluster2): Train-only; Lagged Brand-Cluster Prices
#     - New (k=8), Used (k=7)
# ------------------------------
to_chr <- function(x) as.character(x)

# USED brand clusters (k=7)
used_brand_clust <- fit_clusters_train(df_used_train |> mutate(make = to_chr(make)), by_col = "make", k = 7)
df_used_train <- df_used_train |>
  mutate(make = to_chr(make)) |>
  left_join(used_brand_clust$mapping |> rename(cluster2 = cluster), by = "make")
df_used_test <- df_used_test |>
  mutate(make = to_chr(make)) |>
  left_join(assign_clusters(df_used_test |> mutate(make = to_chr(make)), used_brand_clust) |>
              rename(cluster2 = cluster), by = "make")

used_brand_lag <- lagged_avg_by_group_train(
  df_train   = df_used_train,
  group_cols = "cluster2",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_brand_cluster",
  output_name_lag     = "avg_price_previous_period_brand_cluster"
)

df_used_train <- df_used_train |>
  left_join(used_brand_lag, by = c("yr_month_dtype", "cluster2"))
df_used_test <- df_used_test |>
  left_join(used_brand_lag, by = c("yr_month_dtype", "cluster2"))

# NEW brand clusters (k=8)
new_brand_clust <- fit_clusters_train(df_new_train |> mutate(make = to_chr(make)), by_col = "make", k = 8)
df_new_train <- df_new_train |>
  mutate(make = to_chr(make)) |>
  left_join(new_brand_clust$mapping |> rename(cluster2 = cluster), by = "make")
df_new_test <- df_new_test |>
  mutate(make = to_chr(make)) |>
  left_join(assign_clusters(df_new_test |> mutate(make = to_chr(make)), new_brand_clust) |>
              rename(cluster2 = cluster), by = "make")

new_brand_lag <- lagged_avg_by_group_train(
  df_train   = df_new_train,
  group_cols = "cluster2",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_brand_cluster",
  output_name_lag     = "avg_price_previous_period_brand_cluster"
)

df_new_train <- df_new_train |>
  left_join(new_brand_lag, by = c("yr_month_dtype", "cluster2"))
df_new_test <- df_new_test |>
  left_join(new_brand_lag, by = c("yr_month_dtype", "cluster2"))

# ------------------------------
# 8) Make-Model lag price (train-only)
# ------------------------------
used_mm_lag <- lagged_avg_by_group_train(
  df_train   = df_used_train,
  group_cols = "make_model",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_mm",
  output_name_lag     = "avg_price_previous_period_mm"
)
df_used_train <- df_used_train |> left_join(used_mm_lag, by = c("yr_month_dtype", "make_model"))
df_used_test  <- df_used_test  |> left_join(used_mm_lag, by = c("yr_month_dtype", "make_model"))

new_mm_lag <- lagged_avg_by_group_train(
  df_train   = df_new_train,
  group_cols = "make_model",
  price_col  = "prices_deflated",
  date_col   = "yr_month_dtype",
  lag_months = 12,
  output_name_current = "avg_price_current_period_mm",
  output_name_lag     = "avg_price_previous_period_mm"
)
df_new_train <- df_new_train |> left_join(new_mm_lag, by = c("yr_month_dtype", "make_model"))
df_new_test  <- df_new_test  |> left_join(new_mm_lag, by = c("yr_month_dtype", "make_model"))

# ------------------------------
# 9) PCA (Train-only), then project to Test
# ------------------------------
spec_cols <- c("Engine_Displacement_Cleaned","Engine_Torque_Cleaned","Engine_Capacity","Top_Speed",
               "Acceleration","Dimensions_Length","Dimensions_Width","Dimensions_Wheelbase",
               "Cargo_Volume","Unladen_Weight","Fuel_Economy")

# NEW PCA
X_new_tr  <- scale(df_new_train |> dplyr::select(all_of(spec_cols)))
pca_new   <- princomp(X_new_tr)

df_new_train <- df_new_train |>
  mutate(pc1 = pca_new$scores[,1], pc2 = pca_new$scores[,2])

mu_new <- attr(X_new_tr, "scaled:center")
sd_new <- attr(X_new_tr, "scaled:scale"); sd_new[sd_new==0] <- 1

X_new_te <- scale(df_new_test |> dplyr::select(all_of(spec_cols)), center = mu_new, scale = sd_new)
scores_new_te <- predict(pca_new, newdata = as.data.frame(X_new_te))
df_new_test <- df_new_test |>
  mutate(pc1 = scores_new_te[,1], pc2 = scores_new_te[,2])

# USED PCA
X_used_tr <- scale(df_used_train |> dplyr::select(all_of(spec_cols)))
pca_used  <- princomp(X_used_tr)

df_used_train <- df_used_train |>
  mutate(pc1 = pca_used$scores[,1], pc2 = pca_used$scores[,2])

mu_used <- attr(X_used_tr, "scaled:center")
sd_used <- attr(X_used_tr, "scaled:scale"); sd_used[sd_used==0] <- 1

X_used_te <- scale(df_used_test |> dplyr::select(all_of(spec_cols)), center = mu_used, scale = sd_used)
scores_used_te <- predict(pca_used, newdata = as.data.frame(X_used_te))
df_used_test <- df_used_test |>
  mutate(pc1 = scores_used_te[,1], pc2 = scores_used_te[,2])

# ------------------------------
# 10) Models
# ------------------------------

# -- Baselines (no clusters), using make-model lag (mm)
new_model <- lm(
  prices_deflated ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed +
    Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy +
    Body_type + province + make + mnth + age + yr_decimal,
  data = df_new_train
)
used_model <- lm(
  log(prices_deflated) ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed +
    Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy +
    Body_type + province + make + mnth + age + yr_decimal,
  data = df_used_train
)

# -- PCA variants (no clusters)
new_model_pca <- lm(
  prices_deflated ~ pc1 + pc2 + Body_type + province + make + mnth + age + yr_decimal,
  data = df_new_train
)
used_model_pca <- lm(
  log(prices_deflated) ~ pc1 + pc2 + Body_type + province + make + mnth + age + yr_decimal,
  data = df_used_train
)

# -- Cluster models (make-model clusters), USE ONLY LAGGED CLUSTER PRICE (no current-period feature)
new_model_cluster <- lm(
  prices_deflated ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed +
    Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy +
    Body_type + province + cluster + mnth + age + yr_decimal,
  data = df_new_train
)
used_model_cluster <- lm(
  log(prices_deflated) ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed +
    Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy +
    Body_type + province + cluster + mnth + age + yr_decimal,
  data = df_used_train
)

# -- Cluster + PCA (brand clusters cluster2)
new_model_cluster_pca <- lm(
  prices_deflated ~ pc1 + pc2 + Body_type + province + cluster2 + mnth + age + yr_decimal,
  data = df_new_train
)
used_model_cluster_pca <- lm(
  log(prices_deflated) ~ pc1 + pc2 + Body_type + province + cluster2 + mnth + age + yr_decimal,
  data = df_used_train
)


# Reference Models
new_model_ref <- lm(prices_deflated ~ Body_type + province + make + mnth + age + yr_decimal, data = df_new_train)
used_model_ref <- lm(log(prices_deflated) ~ Body_type + province + make + mnth + age + yr_decimal, data = df_used_train)



#### No fixed effects
new_model_ref2 <- lm(prices_deflated ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed + 
                       Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy, data = df_new_train)
used_model_ref2 <- lm(log(prices_deflated) ~ Engine_Displacement_Cleaned + Engine_Torque_Cleaned + Engine_Capacity + Top_Speed + 
                        Acceleration + Dimensions_Length + Dimensions_Width + Dimensions_Wheelbase + Cargo_Volume + Unladen_Weight + Fuel_Economy, data = df_used_train)


# ------------------------------
# 11) Evaluation
# ------------------------------
rmse_new  <- rmse(predict(new_model, newdata = df_new_test), df_new_test$prices_deflated)
rmse_used <- rmse(exp(predict(used_model, newdata = df_used_test)), df_used_test$prices_deflated)

rmse_new_pca  <- rmse(predict(new_model_pca, newdata = df_new_test), df_new_test$prices_deflated)
rmse_used_pca <- rmse(exp(predict(used_model_pca, newdata = df_used_test)), df_used_test$prices_deflated)

rmse_new_cluster  <- rmse(predict(new_model_cluster, newdata = df_new_test), df_new_test$prices_deflated)
rmse_used_cluster <- rmse(exp(predict(used_model_cluster, newdata = df_used_test)), df_used_test$prices_deflated)

rmse_new_cluster_pca  <- rmse(predict(new_model_cluster_pca, newdata = df_new_test), df_new_test$prices_deflated)
rmse_used_cluster_pca <- rmse(exp(predict(used_model_cluster_pca, newdata = df_used_test)), df_used_test$prices_deflated)

test_predictions_ref_new <- predict(new_model_ref, newdata = df_new_test)
rmse_new_ref <- rmse(test_predictions_ref_new, df_new_test$prices_deflated)
test_predictions_ref_used <- predict(used_model_ref, newdata = df_used_test)
rmse_used_ref <- rmse(exp(test_predictions_ref_used), df_used_test$prices_deflated)

test_predictions_ref2_new <- predict(new_model_ref2, newdata = df_new_test)
rmse_new_ref2 <- rmse(test_predictions_ref2_new, df_new_test$prices_deflated)
test_predictions_ref2_used <- predict(used_model_ref2, newdata = df_used_test)
rmse_used_ref2 <- rmse(exp(test_predictions_ref2_used), df_used_test$prices_deflated)


results_rmse <- tibble(
  model = c("new_model","used_model","new_model_pca","used_model_pca",
            "new_model_cluster","used_model_cluster",
            "new_model_cluster_pca","used_model_cluster_pca", 
            "new_model_ref", "used_model_ref", "new_model_ref2", "used_model_ref2"),
  rmse  = c(rmse_new, rmse_used, rmse_new_pca, rmse_used_pca,
            rmse_new_cluster, rmse_used_cluster,
            rmse_new_cluster_pca, rmse_used_cluster_pca, 
            rmse_new_ref, rmse_used_ref, rmse_new_ref2, rmse_used_ref2)
)

print(results_rmse)

r_squared_tibble <- tibble(
  new_model = summary(new_model)$r.squared,
  new_model_cluster = summary(new_model_cluster)$r.squared,
  # new_model_cluster2 = summary(new_model_cluster2)$r.squared,
  new_model_pca = summary(new_model_pca)$r.squared,
  new_model_cluster_pca = summary(new_model_cluster_pca)$r.squared,
  new_model_ref = summary(new_model_ref)$r.squared,
  new_model_ref2 = summary(new_model_ref2)$r.squared,
  used_model = summary(used_model)$r.squared,
  used_model_cluster = summary(used_model_cluster)$r.squared,
  # used_model_cluster2 = summary(used_model_cluster2)$r.squared,
  used_model_pca = summary(used_model_pca)$r.squared,
  used_model_cluster_pca = summary(used_model_cluster_pca)$r.squared,
  used_model_ref = summary(used_model_ref)$r.squared,
  used_model_ref2 = summary(used_model_ref2)$r.squared
) |>
  pivot_longer(
    cols = everything(),
    names_to = "model",
    values_to = "r_squared"
  )

results <- inner_join(results_rmse, r_squared_tibble)



# Optionally, write outputs
# write_csv(results_rmse, "results_rmse_leakage_safe.csv")
# write_csv(r_squared_tibble, "r2_tibble_leakage_safe.csv")
write_csv2(results, "rmse_r2.csv")
# saveRDS(list(new_model=new_model, used_model=used_model,
#              new_model_pca=new_model_pca, used_model_pca=used_model_pca,
#              new_model_cluster=new_model_cluster, used_model_cluster=used_model_cluster,
#              new_model_cluster_pca=new_model_cluster_pca, used_model_cluster_pca=used_model_cluster_pca),
#        file = "models_leakage_safe.rds")
