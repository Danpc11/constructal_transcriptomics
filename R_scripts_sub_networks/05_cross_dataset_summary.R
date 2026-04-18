#!/usr/bin/env Rscript
# =============================================================================
# 05_cross_dataset_summary.R
# Consolidate metrics and gene-driver scores into cross-dataset summary tables.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# BRANCH: Shared
#
# CORRECTIONS:
#   [C18] trend_summary: physiological ordering via phys_order before
#         first()/last(). Prevents CEI_change inversion due to alphabetical order.
#   [C19] Datasets and state_orders updated to the microarray quartet.
#   [NEW] Warn if the join generates NAs in phys_order (unrecognized conditions).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

# -----------------------------------------------------------------------------
# MAPEO BIOLÓGICO 🔬
# -----------------------------------------------------------------------------
dataset_info <- data.frame(
  accession = c("GSE76895", "GSE18732", "GSE15653", "GSE27951"),
  organ     = c(
    "Pancreatic islets",   # GSE76895: Taneera et al. 2012
    "Skeletal muscle",     # GSE18732: Gallagher et al. 2010
    "Liver",               # GSE15653: Ahrén et al. 2010
    "Adipose tissue"       # GSE27951: Civelek et al. 2017
  ),
  stringsAsFactors = FALSE
)

state_orders <- list(
  GSE76895 = c("ND",  "IGT", "T2D"),
  GSE18732 = c("ND",  "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

message("=== Consolidación con contexto biológico ===")

# -----------------------------------------------------------------------------
# ORDEN FISIOLÓGICO
# -----------------------------------------------------------------------------
order_df <- bind_rows(lapply(targets, function(acc) {
  ord <- state_orders[[acc]]
  data.frame(
    accession  = acc,
    condition  = ord,
    phys_order = seq_along(ord),
    stringsAsFactors = FALSE
  )
}))

# -----------------------------------------------------------------------------
# LECTOR GENÉRICO
# -----------------------------------------------------------------------------
read_and_label <- function(dir_path, pattern, branch_name) {
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) return(NULL)

  bind_rows(lapply(files, function(f) {
    d <- fread(f)
    acc_name <- sub("_.*", "", basename(f))
    if (!"accession" %in% colnames(d)) d$accession <- acc_name
    d$branch <- branch_name
    d
  }))
}

# -----------------------------------------------------------------------------
# 1. MÉTRICAS
# -----------------------------------------------------------------------------
message("1. Métricas...")

metrics <- bind_rows(
  read_and_label("results/metrics", "_constructal_metrics\\.tsv$", "shared"),
  read_and_label("results/full_metrics", "_full_constructal_metrics\\.tsv$", "full")
)

if (nrow(metrics) > 0) {

  metrics_ordered <- metrics |>
    left_join(order_df, by = c("accession", "condition")) |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, phys_order)

  # [C18] Verificar que el join no produjo NAs en phys_order
  # (indicaria condiciones no reconocidas en state_orders)
  na_rows <- metrics_ordered[is.na(metrics_ordered$phys_order), ]
  if (nrow(na_rows) > 0L) {
    warning(nrow(na_rows), " filas con phys_order=NA (condicion no reconocida): ",
            paste(unique(paste0(na_rows$accession, "/", na_rows$condition)),
                  collapse = ", "),
            "\nEstas filas se EXCLUYEN de trend_summary para evitar first()/last() incorrectos.")
    metrics_ordered <- metrics_ordered[!is.na(metrics_ordered$phys_order), ]
  }

  trend_summary <- metrics_ordered |>
    group_by(branch, accession, organ) |>
    summarise(
      healthiest  = first(condition),
      final_state = last(condition),
      EG_change   = last(EG)   - first(EG),
      Hb_change   = last(Hb)   - first(Hb),
      Rbar_change = last(Rbar) - first(Rbar),
      Gbar_change = last(Gbar) - first(Gbar),
      .groups = "drop"
    )

  fwrite(metrics_ordered |> select(-phys_order),
         "results/summary/all_metrics_combined.tsv", sep = "\t")

  fwrite(trend_summary,
         "results/summary/cross_dataset_trends.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 2. DRIVERS
# -----------------------------------------------------------------------------
message("2. Drivers...")

driver_df <- bind_rows(
  read_and_label("results/gene_drivers", "_gene_driver_scores\\.tsv$", "shared"),
  read_and_label("results/full_gene_drivers", "_full_gene_driver_scores\\.tsv$", "full")
)

if (nrow(driver_df) > 0) {

  driver_df <- driver_df |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")"))

  fwrite(driver_df,
         "results/summary/all_gene_driver_scores_combined.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 3. OPTIMUM
# -----------------------------------------------------------------------------
message("3. Optimum...")

opt_df <- bind_rows(
  read_and_label("results/optimum", "_deviation_from_optimum\\.tsv$", "shared"),
  # [CORRECCIÓN] Carpeta cambiada a "results/full_optimum"
  read_and_label("results/full_optimum", "_full_deviation_from_optimum\\.tsv$", "full") 
)

if (nrow(opt_df) > 0) {

  opt_df <- opt_df |>
    left_join(order_df, by = c("accession", "condition")) |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, phys_order) |>
    select(-phys_order)

  fwrite(opt_df,
         "results/summary/all_deviation_from_optimum.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 4. DE INTEGRADO
# -----------------------------------------------------------------------------
message("4. DE integrado...")

de_int_df <- bind_rows(
  # [CORRECCIÓN] Regex anclado al inicio para evitar que cargue los archivos "full"
  read_and_label("results/differential_expression", "^[A-Z0-9]+_integrated_drivers_DE\\.tsv$", "shared"),
  read_and_label("results/differential_expression", "_full_integrated_drivers_DE\\.tsv$", "full")
)

if (nrow(de_int_df) > 0) {

  de_int_df <- de_int_df |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, desc(combined_score))

  fwrite(de_int_df,
         "results/summary/all_integrated_drivers_DE.tsv", sep = "\t")
}

message("=== DONE ===")