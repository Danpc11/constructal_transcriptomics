#!/usr/bin/env Rscript
# =============================================================================
# 02_networks_modularity_constructal_metrics.R
# Construye redes de co-expresion por estado sobre el espacio COMUN de genes
# (interseccion de los 4 datasets de microarray), calcula metricas constructales
# y detecta modulos.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# RAMA: Compartida (genes fijos, metricas directamente comparables entre datasets)
#
# CORRECCIONES:
#   [C5]  Subred fija de genes comunes filtrada por varianza global (script 01).
#   [C6]  CEI calculado globalmente sobre todos los datasets/estados juntos.
#   [C7]  find_modules maximiza n_mod * Q (modularidad Newman) en lugar de solo
#         n_mod. Fallback a corte mediano si ningun corte produce modulos validos.
#   [C8]  Orden fisiologico en archivos de salida.
#   [FIX-Rbar]  Rbar = 2*tr(L+)/p. Version anterior usaba /(p-1) -> error x p.
#   [FIX-match] match(states, condition) en lugar del match invertido anterior.
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(expm)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
dir.create("results/metrics",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/modules",  recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

# Canonical physiological order
state_orders <- list(
  GSE76895 = c("ND", "IGT", "T2D"),
  GSE18732 = c("ND", "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

# -----------------------------------------------------------------------------
# Network functions
# -----------------------------------------------------------------------------

compute_adjacency <- function(expr, cor_method = "bicor") {
  expr_t  <- t(expr)
  sft     <- pickSoftThreshold(expr_t, dataIsExpr = TRUE, corFnc = cor_method,
                               powerVector = 1:20, verbose = 0)
  beta    <- sft$powerEstimate
  if (is.na(beta)) {
    fit       <- sft$fitIndices
    candidate <- fit$Power[fit$SFT.R.sq >= 0.80]
    beta      <- if (length(candidate) > 0) min(candidate) else 6L
  }
  cor_mat <- if (cor_method == "bicor") {
    bicor(expr_t, maxPOutliers = 0.1)
  } else {
    cor(expr_t, method = "pearson")
  }
  cor_mat[is.na(cor_mat)] <- 0
  adj             <- abs(cor_mat)^beta
  diag(adj)       <- 0
  rownames(adj)   <- colnames(adj) <- rownames(expr)
  list(adj = adj, beta = beta, sft = sft)
}

# [C7] Seleccion de corte maximizando n_mod * Q (no solo n_mod)
find_modules <- function(W, min_size = 20L) {
  d  <- as.dist(1 - W)
  hc <- hclust(d, method = "average")
  g  <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                    diag = FALSE)
  best <- NULL
  for (q in seq(0.55, 0.80, by = 0.05)) {
    h       <- quantile(hc$height, probs = q, na.rm = TRUE)
    cl      <- cutree(hc, h = h)
    tab     <- table(cl)
    keep_cl <- names(tab)[tab >= min_size]
    n_mod   <- length(keep_cl)
    if (n_mod < 1L) next
    module_tmp        <- paste0("M", cl)
    small             <- names(table(module_tmp))[table(module_tmp) < min_size]
    module_tmp[module_tmp %in% small] <- "grey"
    mem     <- as.integer(factor(module_tmp))
    Q_val   <- tryCatch(
      modularity(g, membership = mem, weights = E(g)$weight),
      error = function(e) 0
    )
    score <- n_mod * max(Q_val, 0)
    candidate <- list(h = h, cl = cl, n_mod = n_mod, Q = Q_val, score = score)
    if (is.null(best) || candidate$score > best$score) best <- candidate
  }
  if (is.null(best)) {
    # Fallback: corte mediano si ningun percentil produce modulos de tamano >= min_size
    h    <- median(hc$height)
    cl   <- cutree(hc, h = h)
    best <- list(cl = cl, n_mod = 1L, Q = 0, score = 0)
  }
  module        <- paste0("M", best$cl)
  tab           <- table(module)
  small         <- names(tab)[tab < min_size]
  module[module %in% small] <- "grey"
  names(module) <- rownames(W)
  list(module    = module,
       n_modules = length(setdiff(unique(module), "grey")),
       Q         = best$Q)
}

# -----------------------------------------------------------------------------
# Metricas constructales
# [FIX-Rbar] Rbar = 2*tr(L+)/p  (antes /(p-1) -> sobreestimacion por factor p)
# -----------------------------------------------------------------------------

compute_strength_entropy <- function(W) {
  k <- rowSums(W)
  p <- k / sum(k)
  -sum(p * log(p + 1e-12))
}

compute_global_efficiency <- function(W) {
  g    <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
  D    <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / D
  diag(invD) <- NA_real_
  mean(invD[upper.tri(invD)], na.rm = TRUE) * 2
}

compute_communicability <- function(W) {
  k    <- rowSums(W)
  Dinv <- diag(1 / sqrt(k + 1e-12), nrow(W))
  G    <- expm(Dinv %*% W %*% Dinv)
  mean(G[row(G) != col(G)])
}

compute_avg_effective_resistance <- function(W) {
  p        <- nrow(W)
  k        <- rowSums(W)
  L        <- diag(k, p) - W
  eig      <- eigen(L, symmetric = TRUE)
  tol      <- max(abs(eig$values)) * p * 1e-10
  inv_vals <- ifelse(eig$values > tol, 1 / eig$values, 0)
  Lplus    <- eig$vectors %*% diag(inv_vals, p) %*% t(eig$vectors)
  as.numeric(2 * sum(diag(Lplus)) / p)   # [FIX-Rbar]
}

# Wrapper unico — mismo en scripts 03 y 04 para consistencia
compute_metrics <- function(W) {
  c(EG   = compute_global_efficiency(W),
    Gbar = compute_communicability(W),
    Rbar = compute_avg_effective_resistance(W),
    Hb   = compute_strength_entropy(W))
}

# -----------------------------------------------------------------------------
# [C5] Carga del conjunto fijo de genes comunes (generado por script 01)
# -----------------------------------------------------------------------------

hv_genes_path <- "data/processed/high_variance_genes.rds"
if (!file.exists(hv_genes_path)) {
  stop("'high_variance_genes.rds' no encontrado. Ejecuta primero el script 01.")
}
hv_genes    <- readRDS(hv_genes_path)
p_sub       <- min(800L, length(hv_genes))
fixed_genes <- hv_genes[seq_len(p_sub)]
message("Fixed subnetwork nodes: ", length(fixed_genes))

# -----------------------------------------------------------------------------
# Bucle principal
# -----------------------------------------------------------------------------

metrics_list <- list()

for (acc in targets) {
  obj   <- readRDS(file.path("data/processed", paste0(acc, "_processed.rds")))
  expr  <- obj$expr
  pheno <- obj$pheno

  ord    <- state_orders[[acc]]
  states <- ord[ord %in% as.character(unique(pheno$condition))]

  acc_metrics <- list()

  for (st in states) {
    samples  <- pheno$.sample_id[pheno$condition == st]
    expr_sub <- expr[, samples, drop = FALSE]
    if (ncol(expr_sub) < 4L) {
      message("  Skipping ", acc, "/", st, ": too few samples")
      next
    }

    net <- compute_adjacency(expr_sub, cor_method = "bicor")
    W   <- net$adj

    # Subred fija: mismos nodos en todos los datasets/estados  [C5]
    common_avail <- intersect(fixed_genes, rownames(W))
    if (length(common_avail) < 50L) {
      message("  Warning: only ", length(common_avail),
              " fixed genes available in ", acc, "/", st)
    }
    Wsub <- W[common_avail, common_avail, drop = FALSE]

    # Modulos con seleccion por n*Q  [C7]
    mods      <- find_modules(Wsub)
    module_df <- data.frame(gene      = names(mods$module),
                            module    = unname(mods$module),
                            accession = acc,
                            condition = st)
    fwrite(module_df,
           file.path("results/modules",
                     paste0(acc, "_", st, "_modules.tsv")), sep = "\t")

    m <- compute_metrics(Wsub)

    acc_metrics[[st]] <- data.frame(
      accession    = acc,
      condition    = st,
      n_samples    = ncol(expr_sub),
      n_genes      = nrow(expr_sub),
      n_subgenes   = nrow(Wsub),
      beta         = net$beta,
      EG           = m["EG"],
      Gbar         = m["Gbar"],
      Rbar         = m["Rbar"],
      Hb           = m["Hb"],
      mean_k       = mean(rowSums(Wsub)),
      n_modules    = mods$n_modules,
      Q_modularity = mods$Q
    )

    saveRDS(list(W = Wsub, genes = rownames(Wsub), beta = net$beta),
            file.path("results/networks",
                      paste0(acc, "_", st, "_network.rds")))
  }

  # [C8 + FIX-match] Orden fisiologico garantizado.
  # match(states, condition) devuelve el indice de fila de cada estado en el
  # orden correcto. La version incorrecta match(condition, states) devolvia la
  # posicion en el vector states para cada fila — mezclaba en lugar de ordenar.
  metrics_df <- dplyr::bind_rows(acc_metrics)
  row_order  <- match(states, metrics_df$condition)
  metrics_df <- metrics_df[row_order[!is.na(row_order)], ]

  fwrite(metrics_df,
         file.path("results/metrics",
                   paste0(acc, "_constructal_metrics.tsv")), sep = "\t")
  metrics_list[[acc]] <- metrics_df
}

# [C6] CEI calculado GLOBALMENTE sobre todos los datasets/estados
all_metrics <- dplyr::bind_rows(metrics_list)

zscore_global <- function(x) as.numeric(scale(x))

all_metrics <- all_metrics |>
  dplyr::mutate(
    CEI = zscore_global(EG) + zscore_global(Gbar) -
          zscore_global(Rbar) + zscore_global(Hb)
  )

# Reescribir archivos por dataset con CEI global incluido
for (acc in targets) {
  df_acc <- all_metrics[all_metrics$accession == acc, ]
  fwrite(df_acc,
         file.path("results/metrics",
                   paste0(acc, "_constructal_metrics.tsv")), sep = "\t")
}

fwrite(all_metrics, "results/metrics/all_constructal_metrics.tsv", sep = "\t")
message("Done.")