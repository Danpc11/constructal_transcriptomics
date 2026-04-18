#!/usr/bin/env Rscript
# =============================================================================
# 04c_expression_differential.R
# Differential expression with limma for the four microarray datasets and
# integration with gene-driver scores (PTI / IRI / TRI).
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
#
# NOTE: limma is the correct method for microarrays. Data are already
#       RMA-normalized on the log2 scale from script 01. Do NOT use
#       DESeq2 or edgeR (designed for raw RNA-seq counts).
#
# CORRECTIONS vs. the original user version:
#   - Datasets updated to the four microarray datasets of the metabolic quartet.
#     (removed GSE50244 RNA-seq and GSE152991, which are not part of this pipeline)
#   - state_orders identical to the rest of the pipeline (scripts 02-05).
#   - healthy_state defined explicitly per dataset for correct releveling
#     in the linear model (healthy state = reference factor level).
#   - Additional critical-transition contrast: T2D_vs_IGT and ObT2D_vs_ObnoT2D
#     to capture the phase-transition signature at the critical step.
#   - Direct integration with gene_driver_scores.tsv (shared branch) and
#     full_gene_driver_scores.tsv (full branch), generating:
#       combined_score = (PTI + IRI + TRI) * |logFC|
#
# OUTPUTS:
#   results/differential_expression/{acc}_limma_DE.tsv
#   results/differential_expression/all_datasets_limma_DE.tsv
#   results/differential_expression/{acc}_integrated_drivers_DE.tsv
#   results/differential_expression/{acc}_full_integrated_drivers_DE.tsv
# =============================================================================
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

# Orden fisiológico
state_orders <- list(
  GSE76895 = c("ND",  "IGT", "T2D"),
  GSE18732 = c("ND",  "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

# Estado sano de referencia
healthy_state <- list(
  GSE76895 = "ND",
  GSE18732 = "ND",
  GSE15653 = "Lean",
  GSE27951 = "NGT"
)

# -----------------------------------------------------------------------------
# Funciones limma robustecidas
# -----------------------------------------------------------------------------

run_limma <- function(expr, pheno, states, ref_state, acc) {
  
  # Filtrar y asegurar alineación
  pheno <- dplyr::filter(pheno, condition %in% states)
  expr  <- expr[, pheno$.sample_id, drop = FALSE]

  # [CORREGIDO] Factorizar asegurando que el ref_state sea el nivel base
  pheno$condition <- factor(pheno$condition, levels = states)
  pheno$condition <- relevel(pheno$condition, ref = ref_state)

  # Matriz de diseño (sin intercepto para facilitar contrastes)
  design <- model.matrix(~ 0 + condition, data = pheno)
  
  # [CORREGIDO] Eliminar el prefijo "condition" de forma segura
  colnames(design) <- gsub("condition", "", colnames(design))

  fit <- lmFit(expr, design)

  # --- Generación dinámica de contrastes ---
  contrasts_list <- list()

  # 1. Todos contra el estado sano
  non_ref <- setdiff(colnames(design), ref_state)
  for (st in non_ref) {
    cname <- paste0(st, "_vs_", ref_state)
    contrasts_list[[cname]] <- paste0(st, " - ", ref_state)
  }

  # 2. Transiciones críticas (Intermedio -> Final)
  if (all(c("IGT", "T2D") %in% colnames(design)))
    contrasts_list[["T2D_vs_IGT"]] <- "T2D - IGT"

  if (all(c("Obese_noT2D", "Obese_T2D") %in% colnames(design)))
    contrasts_list[["ObT2D_vs_ObnoT2D"]] <- "Obese_T2D - Obese_noT2D"

  if (length(contrasts_list) == 0L) {
    message("  -> No hay contrastes válidos para ", acc)
    return(NULL)
  }

  contrast_matrix <- makeContrasts(contrasts = unlist(contrasts_list), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

  # Extraer resultados para todos los contrastes
  results <- lapply(colnames(contrast_matrix), function(coef_name) {
    tt           <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
    tt$gene      <- rownames(tt)
    tt$contrast  <- coef_name
    tt$accession <- acc
    tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "contrast", "accession")]
  })

  dplyr::bind_rows(results)
}

# -----------------------------------------------------------------------------
# Bucle Principal de Expresión Diferencial
# -----------------------------------------------------------------------------

all_results <- list()

for (acc in targets) {
  message("======================================")
  message("DE limma: ", acc)

  proc_file <- file.path("data/processed", paste0(acc, "_processed_full.rds"))
  if (!file.exists(proc_file)) {
    proc_file <- file.path("data/processed", paste0(acc, "_processed.rds"))
  }
  if (!file.exists(proc_file)) {
    message("  -> Archivo procesado no encontrado, saltando.")
    next
  }

  obj   <- readRDS(proc_file)
  expr  <- obj$expr
  pheno <- obj$pheno

  # Obtener los estados que realmente existen en el fenotipo
  states    <- state_orders[[acc]]
  states    <- states[states %in% as.character(unique(pheno$condition))]
  ref_state <- healthy_state[[acc]]

  if (length(states) < 2L || !ref_state %in% states) {
    message("  -> Estados insuficientes o falta referencia. Saltando.")
    next
  }

  res <- tryCatch({
    run_limma(expr, pheno, states, ref_state, acc)
  }, error = function(e) {
    message("  -> Error en limma: ", e$message)
    NULL
  })

  if (!is.null(res)) {
    out_path <- file.path("results/differential_expression", paste0(acc, "_limma_DE.tsv"))
    fwrite(res, out_path, sep = "\t")
    message("  -> ", nrow(res), " filas | ", length(unique(res$contrast)), 
            " contrastes | guardado: ", basename(out_path))
    all_results[[acc]] <- res
  }
}

if (length(all_results) > 0L) {
  all_df <- dplyr::bind_rows(all_results)
  fwrite(all_df, "results/differential_expression/all_datasets_limma_DE.tsv", sep = "\t")
  message("Consolidado global DE guardado.")
}

# -----------------------------------------------------------------------------
# Integración de Scores (Topología + DE)
# -----------------------------------------------------------------------------

integrate_drivers <- function(de_file, driver_file, out_file) {
  if (!file.exists(de_file) || !file.exists(driver_file)) {
    message("  -> Faltan archivos para integración: ", basename(de_file), " o ", basename(driver_file))
    return(invisible(NULL))
  }

  de <- fread(de_file)
  dr <- fread(driver_file)

  # Extraer el mejor contraste (mayor |logFC|) por gen
  de_best <- de |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(gene, logFC, adj.P.Val, contrast)

  driver_cols <- intersect(c("gene", "PTI", "IRI", "TRI", "KO_support", "class"), colnames(dr))

  # Normalización rank-based (percentil) para combined_score.
  # min_max_scale tiene un artefacto: el gen con score minimo siempre obtiene
  # combined_score=0 sin importar su logFC. rank_scale lo evita:
  # convierte scores a percentiles [0,1] preservando el orden relativo
  # sin distorsionar los extremos.
  rank_scale <- function(x) {
    r <- rank(x, na.last = "keep", ties.method = "average")
    r / max(r, na.rm = TRUE)
  }

  integrated <- dr |>
    dplyr::select(dplyr::all_of(driver_cols)) |>
    dplyr::inner_join(de_best, by = "gene") |>
    dplyr::mutate(
      DE_sig = adj.P.Val < 0.05,
      # Suma topológica — sustituir NA por 0 si alguna columna falta
      PTI = if ("PTI" %in% colnames(dr)) PTI else 0,
      IRI = if ("IRI" %in% colnames(dr)) IRI else 0,
      TRI = if ("TRI" %in% colnames(dr)) TRI else 0,
      raw_driver_sum = rowSums(cbind(PTI, IRI, TRI), na.rm = TRUE),
      # Percentil en [0,1]: preserva orden sin artefacto de borde
      scaled_driver  = rank_scale(raw_driver_sum),
      # Score final: impacto topologico (percentil) x magnitud DE
      combined_score = scaled_driver * abs(logFC)
    ) |>
    dplyr::arrange(dplyr::desc(combined_score))

  fwrite(integrated, out_file, sep = "\t")

  n_tr <- sum(integrated$class == "transition_rescue" & integrated$DE_sig, na.rm = TRUE)
  message("  -> ", basename(out_file), " | ", nrow(integrated), " genes | ",
          n_tr, " transition_rescue significativos (adj.P < 0.05)")

  invisible(integrated)
}

message("======================================")
message("Integrando DE con gene drivers (Rama Shared)...")
for (acc in targets) {
  integrate_drivers(
    de_file     = file.path("results/differential_expression", paste0(acc, "_limma_DE.tsv")),
    driver_file = file.path("results/gene_drivers", paste0(acc, "_gene_driver_scores.tsv")),
    out_file    = file.path("results/differential_expression", paste0(acc, "_integrated_drivers_DE.tsv"))
  )
}

message("Integrando DE con gene drivers (Rama Full)...")
for (acc in targets) {
  integrate_drivers(
    de_file     = file.path("results/differential_expression", paste0(acc, "_limma_DE.tsv")),
    driver_file = file.path("results/full_gene_drivers", paste0(acc, "_full_gene_driver_scores.tsv")),
    out_file    = file.path("results/differential_expression", paste0(acc, "_full_integrated_drivers_DE.tsv"))
  )
}

message("Pipeline Completado.")