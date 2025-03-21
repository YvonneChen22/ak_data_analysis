common DEG list
================
Sida Chen

``` r
ak1_up <- read.csv("./DEGs/ak1_7_up.csv") |>
  janitor::clean_names() |>
  filter(cluster_4_p_value < 0.05) |>
  dplyr::rename(
    ak1_dys_average = cluster_4_average,
    ak1_dys_log2_fold_change = cluster_4_log2_fold_change,
    ak1_dys_p_value = cluster_4_p_value,
    ak1_norm_average = cluster_6_average,
    ak1_norm_log2_fold_change = cluster_6_log2_fold_change,
    ak1_norm_p_value = cluster_6_p_value
  ) |>
  select(-"feature_id")
  

ak1_down <- read.csv("./DEGs/ak1_7_down.csv") |>
  janitor::clean_names() |>
  filter(cluster_4_p_value < 0.05) |>
  dplyr::rename(
    ak1_dys_average = cluster_4_average,
    ak1_dys_log2_fold_change = cluster_4_log2_fold_change,
    ak1_dys_p_value = cluster_4_p_value,
    ak1_norm_average = cluster_6_average,
    ak1_norm_log2_fold_change = cluster_6_log2_fold_change,
    ak1_norm_p_value = cluster_6_p_value
  ) |>
  select(-"feature_id")

ak2_up <- read.csv("./DEGs/ak2_3_up.csv") |>
  janitor::clean_names() |>
  filter(cluster_1_p_value < 0.05) |>
  dplyr::rename(
    ak2_dys_average = cluster_2_average,
    ak2_dys_log2_fold_change = cluster_2_log2_fold_change,
    ak2_dys_p_value = cluster_2_p_value,
    ak2_norm_average = cluster_1_average,
    ak2_norm_log2_fold_change = cluster_1_log2_fold_change,
    ak2_norm_p_value = cluster_1_p_value
  ) |>
  select(-"feature_id")

ak2_down <- read.csv("./DEGs/ak2_3_down.csv") |>
  janitor::clean_names() |>
  filter(cluster_1_p_value < 0.05) |>
 dplyr::rename(
    ak2_dys_average = cluster_2_average,
    ak2_dys_log2_fold_change = cluster_2_log2_fold_change,
    ak2_dys_p_value = cluster_2_p_value,
    ak2_norm_average = cluster_1_average,
    ak2_norm_log2_fold_change = cluster_1_log2_fold_change,
    ak2_norm_p_value = cluster_1_p_value
  ) |>
  select(-"feature_id")

ak3_up <- read.csv("./DEGs/ak3_5_up.csv") |>
  janitor::clean_names() |>
  filter(cluster_4_p_value < 0.05) |>
  dplyr::rename(
    ak3_dys_average = cluster_4_average,
    ak3_dys_log2_fold_change = cluster_4_log2_fold_change,
    ak3_dys_p_value = cluster_4_p_value,
    ak3_norm_average = cluster_3_average,
    ak3_norm_log2_fold_change = cluster_3_log2_fold_change,
    ak3_norm_p_value = cluster_3_p_value
  ) |>
  select(-"feature_id")

ak3_down <- read.csv("./DEGs/ak3_5_down.csv") |>
  janitor::clean_names() |>
  filter(cluster_4_p_value < 0.05) |>
  dplyr::rename(
    ak3_dys_average = cluster_4_average,
    ak3_dys_log2_fold_change = cluster_4_log2_fold_change,
    ak3_dys_p_value = cluster_4_p_value,
    ak3_norm_average = cluster_3_average,
    ak3_norm_log2_fold_change = cluster_3_log2_fold_change,
    ak3_norm_p_value = cluster_3_p_value
  ) |>
  select(-"feature_id")

ak4_up <- read.csv("./DEGs/ak4_7_up.csv") |>
  janitor::clean_names() |>
  filter(cluster_3_p_value < 0.05) |>
  dplyr::rename(
    ak4_dys_average = cluster_2_average,
    ak4_dys_log2_fold_change = cluster_2_log2_fold_change,
    ak4_dys_p_value = cluster_2_p_value,
    ak4_norm_average = cluster_3_average,
    ak4_norm_log2_fold_change = cluster_3_log2_fold_change,
    ak4_norm_p_value = cluster_3_p_value
  ) |>
  select(-"feature_id")

ak4_down <- read.csv("./DEGs/ak4_7_down.csv") |>
  janitor::clean_names() |>
  filter(cluster_3_p_value < 0.05) |>
  dplyr::rename(
    ak4_dys_average = cluster_2_average,
    ak4_dys_log2_fold_change = cluster_2_log2_fold_change,
    ak4_dys_p_value = cluster_2_p_value,
    ak4_norm_average = cluster_3_average,
    ak4_norm_log2_fold_change = cluster_3_log2_fold_change,
    ak4_norm_p_value = cluster_3_p_value
  ) |>
  select(-"feature_id")
```

# Common DEG lists

``` r
# Up-regulated
up_common <- ak1_up |>
  inner_join(ak2_up, by = "feature_name") |>
  inner_join(ak3_up, by = "feature_name") |>
  inner_join(ak4_up, by = "feature_name")

# Down-regulated
down_common <- ak1_down |>
  inner_join(ak2_down, by = "feature_name") |>
  inner_join(ak3_down, by = "feature_name") |>
  inner_join(ak4_down, by = "feature_name")
```

``` r
write_xlsx(up_common, "up_common.xlsx")
write_xlsx(down_common, "down_common.xlsx")
```
