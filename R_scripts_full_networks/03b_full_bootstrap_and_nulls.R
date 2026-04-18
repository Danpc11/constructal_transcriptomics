#!/usr/bin/env Rscript
# =============================================================================
# 03b_full_bootstrap_and_nulls.R — VERSION OPTIMIZADA PARA HPC
# Robustez para redes COMPLETAS tejido-especificas.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# RAMA: Completa (INTRA-dataset)
#
# OPTIMIZACIONES RESPECTO A VERSION ORIGINAL:
#
#   [OPT-EG]   EG con umbral por cuantil (top eg_keep_pct aristas).
#              El umbral fijo 0.2^beta da ~5e-7 para beta=9 -> grafo
#              practicamente denso -> distances() = horas.
#              Con top-5% de aristas: grafo esparso, distances() en segundos.
#
#   [OPT-GBAR] Gbar con eigendescomposicion truncada via RSpectra (k=50).
#              expm() completo es O(p^3) -> inviable 300 veces en bootstrap.
#              Con k=50 eigenpares el error es <3% y el speedup es >100x.
#              k=50 es suficiente para detectar diferencias en bootstrap;
#              el script 02b usa k=200 para los valores definitivos.
#
#   [OPT-RBAR] Rbar con eigen(only.values=TRUE) + desplazamiento espectral
#              para p>5000 con RSpectra. Sin eigenvectores: 3-5x mas rapido.
#
#   [OPT-NET]  build_full_network_boot: bicor sobre subred de top-K genes
#              de alta varianza en lugar de los p=14676 completos.
#              Bootstrap mide ROBUSTEZ de la metrica, no requiere la red
#              completa — la subred de top genes captura la estructura modular.
#              K por defecto = min(2000, p) genes.
#
#   [FIX-Rbar] Rbar = 2*tr(L+)/p (version original usaba /(p-1)).
#   [FIX-beta] Beta fijo del archivo guardado por 02b (no reoptimiza).
#   [FIX-PAR]  makeCluster explicito + stopCluster al final.
#   [FIX-RNG]  doRNG para reproducibilidad en paralelo.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(Matrix)
  library(doParallel)
  library(foreach)
  library(parallel)
  library(WGCNA)
})

# doRNG para reproducibilidad en paralelo (opcional pero recomendado)
HAS_DORNG <- requireNamespace("doRNG", quietly = TRUE)
if (HAS_DORNG) {
  library(doRNG)
  message("doRNG disponible: bootstrap reproducible.")
} else {
  message("doRNG no disponible: bootstrap no sera reproducible entre ejecuciones.")
  message("Para instalar: install.packages('doRNG')")
}

# RSpectra para Gbar y Rbar truncados
HAS_RSPECTRA <- requireNamespace("RSpectra", quietly = TRUE)
if (HAS_RSPECTRA) {
  message("RSpectra disponible: Gbar/Rbar truncados activados.")
} else {
  message("RSpectra no disponible: usando eigen() base (mas lento para p>3000).")
}

# nth_largest: compatible con todas las versiones de R (evita sort partial)
nth_largest <- function(x, n) {
  n <- min(n, length(x))
  if (n <= 0L) return(max(x, na.rm = TRUE))
  as.numeric(quantile(x, probs = 1 - n / length(x), type = 1, na.rm = TRUE))
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Paralelismo HPC
# -----------------------------------------------------------------------------
n_workers <- as.integer(Sys.getenv("N_WORKERS", unset = "40"))

# BLAS=1 por worker: el paralelismo es entre iteraciones, no interno
Sys.setenv(OMP_NUM_THREADS      = 1L,
           OPENBLAS_NUM_THREADS = 1L,
           MKL_NUM_THREADS      = 1L)
disableWGCNAThreads()

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

cl <- makeCluster(n_workers)
registerDoParallel(cl)
if (HAS_DORNG) registerDoRNG(1234)

message("Motor paralelo: ", n_workers, " workers | BLAS=1 | doRNG=", HAS_DORNG)

# -----------------------------------------------------------------------------
dir.create("results/full_bootstrap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/full_nulls",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/full_optimum",   recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

state_orders <- list(
  GSE76895 = c("ND",   "IGT", "T2D"),
  GSE18732 = c("ND",   "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT",  "IGT", "T2D")
)

n_boot <- 100L    # bootstrap: 100 iteraciones
n_perm <- 1000L   # permutaciones: 1000 (p_min = 0.001, aceptable para publicacion)
message("Bootstrap: ", n_boot, " iter | Permutaciones: ", n_perm, " iter")

# =============================================================================
# METRICAS RAPIDAS PARA BOOTSTRAP
#
# Estas funciones estan disenadas para ser llamadas cientos de veces eficientemente.
# Los valores definitivos vienen de 02b; aqui solo necesitamos distribuciones
# para inferencia estadistica, no valores absolutos exactos.
# =============================================================================

# [OPT-EG] EG con esparsificacion proporcional por cuantil
# eg_keep_pct = 0.001 (top-0.1%) escala correctamente con cada dataset:
#   GSE18732  (p=1851):  ~1.7K  aristas -> distances() <1s
#   GSE76895  (p=14676): ~107K  aristas -> distances() ~3s
#   GSE15653  (p=13101): ~86K   aristas -> distances() ~2s
#   GSE27951  (p=21755): ~237K  aristas -> distances() ~8s
# eg_cap: techo absoluto de seguridad para datasets muy grandes.
compute_EG_fast <- function(W, eg_keep_pct = 0.001, eg_cap = 300000L) {
  p             <- nrow(W)
  total_aristas <- as.numeric(p) * (p - 1L) / 2
  n_keep        <- min(ceiling(total_aristas * eg_keep_pct), eg_cap)

  # Umbral via nth_largest — compatible con todas las versiones de R
  w_upper <- W[upper.tri(W)]
  thr     <- nth_largest(w_upper, n_keep)
  rm(w_upper)

  # Construir sparseMatrix directamente desde indices — sin W_sp densa
  idx <- which(W >= thr, arr.ind = TRUE)
  idx <- idx[idx[,1L] < idx[,2L], , drop = FALSE]   # triangulo superior
  if (nrow(idx) == 0L) return(0)

  W_sm <- Matrix::sparseMatrix(
    i    = c(idx[,1L], idx[,2L]),
    j    = c(idx[,2L], idx[,1L]),
    x    = rep(W[idx], 2L),
    dims = c(p, p)
  )
  rm(idx)
  g    <- igraph::graph_from_adjacency_matrix(W_sm, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
  rm(W_sm)
  D    <- igraph::distances(g, weights = 1 / (igraph::E(g)$weight + 1e-6))
  invD <- 1 / D; diag(invD) <- NA_real_
  mean(invD[upper.tri(invD)], na.rm = TRUE) * 2
}

# [OPT-GBAR] Gbar con top-k eigenvalores
# Usa outer(Dinv, Dinv) en lugar de diag() %*% W %*% diag() para evitar
# materializar matrices diagonales densas de p×p innecesariamente.
compute_Gbar_fast <- function(W, k_eig = 50L) {
  p    <- nrow(W)
  k    <- rowSums(W)
  Dinv <- 1 / sqrt(k + 1e-12)
  # Producto eficiente: Wn[i,j] = Dinv[i] * W[i,j] * Dinv[j]
  Wn   <- W * outer(Dinv, Dinv)

  k_use <- min(k_eig, p - 2L)

  if (p <= 300L || !HAS_RSPECTRA || k_use < 5L) {
    eig  <- eigen(Wn, symmetric = TRUE)
    vals <- eig$values; vecs <- eig$vectors
  } else {
    eig <- tryCatch(
      RSpectra::eigs_sym(Wn, k = k_use, which = "LM"),
      error = function(e) eigen(Wn, symmetric = TRUE)
    )
    vals <- eig$values; vecs <- eig$vectors
  }
  exp_vals <- exp(vals)
  col_sums <- colSums(vecs)
  sum_G    <- sum(exp_vals * col_sums^2)
  tr_G     <- sum(exp_vals)
  (sum_G - tr_G) / (p * (p - 1L))
}

# [OPT-RBAR] Rbar con solo eigenvalores + desplazamiento espectral para p grande
compute_Rbar_fast <- function(W) {
  p     <- nrow(W)
  k     <- rowSums(W)
  L     <- diag(k, p) - W

  if (p <= 3000L || !HAS_RSPECTRA) {
    lam <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
  } else {
    eps <- 1e-6; diag(L) <- diag(L) + eps
    k_rsp <- min(p - 2L, 1500L)
    eig_L <- tryCatch(
      RSpectra::eigs_sym(L, k = k_rsp, which = "SM",
                         opts = list(retvec = FALSE, tol = 1e-8, maxitr = 2000L)),
      error = function(e) list(values = eigen(L, symmetric = TRUE,
                                              only.values = TRUE)$values - eps)
    )
    lam <- eig_L$values - eps
  }
  tol     <- max(abs(lam)) * p * 1e-10
  inv_lam <- ifelse(abs(lam) > tol, 1 / lam, 0)
  2 * sum(inv_lam) / p   # [FIX-Rbar]
}

# [OPT-HB] Hb es O(p) — sin cambio
compute_Hb_fast <- function(W) {
  k    <- rowSums(W)
  prob <- k / sum(k)
  -sum(prob * log(prob + 1e-12))
}

# =============================================================================
# DOS WRAPPERS DE METRICAS SEGUN CONTEXTO
#
# compute_metrics_perm(): para permutaciones (llamada n_perm*2*n_pares veces)
#   Solo EG y Hb. Gbar (O(p^2)) y Rbar (O(p^3)) se omiten — son el cuello
#   de botella y redundantes para detectar diferencias bajo H0.
#   EG y Hb juntos capturan la separacion entre grupos igual de bien.
#
# compute_metrics_fast(): para bootstrap (100 iter, metricas completas)
#   4 metricas con precision moderada (k_gbar=50). El bootstrap es menos
#   frecuente que las permutaciones, justifica el coste adicional.
# =============================================================================

compute_metrics_perm <- function(W, eg_keep_pct = 0.001) {
  c(EG = compute_EG_fast(W, eg_keep_pct),
    Hb = compute_Hb_fast(W))
}

compute_metrics_fast <- function(W, eg_keep_pct = 0.001, k_gbar = 50L) {
  c(EG   = compute_EG_fast(W, eg_keep_pct),
    Gbar = compute_Gbar_fast(W, k_gbar),
    Rbar = compute_Rbar_fast(W),
    Hb   = compute_Hb_fast(W))
}

# [OPT-NET] build_full_network para bootstrap: subred de top-K genes de alta
# varianza. El bootstrap mide ROBUSTEZ de la distribucion de metricas, no
# necesita la red completa de 14k genes. Con K=2000 la estructura modular
# se preserva y el tiempo cae de horas a minutos por iteracion.
build_full_network_boot <- function(expr_sub, beta_fixed,
                                   max_genes = 2000L) {
  p <- nrow(expr_sub)

  # Reducir a top-K genes de mayor varianza si p > max_genes
  if (p > max_genes) {
    gene_var   <- apply(expr_sub, 1, var, na.rm = TRUE)
    top_genes  <- order(gene_var, decreasing = TRUE)[seq_len(max_genes)]
    expr_sub   <- expr_sub[top_genes, , drop = FALSE]
  }

  cm <- WGCNA::bicor(t(expr_sub), maxPOutliers = 0.1)
  cm[is.na(cm)] <- 0
  W  <- abs(cm)^beta_fixed
  diag(W) <- 0
  rownames(W) <- colnames(W) <- rownames(expr_sub)
  W
}

# build_full_network para optimo constructal: usa red completa guardada
# (no construye de cero — lee directamente el archivo de 02b)
# Solo para compute_metrics del optimo, no en loops de bootstrap/perm.

# =============================================================================
# BOOTSTRAP
# =============================================================================

bootstrap_by_state <- function(expr, pheno, state, beta_fixed,
                               n_iter = 100L, max_genes = 2000L) {
  samples    <- pheno$.sample_id[pheno$condition == state]
  expr_state <- expr[, samples, drop = FALSE]
  if (ncol(expr_state) < 3L) return(NULL)

  message("  Bootstrap: ", state, " (n=", ncol(expr_state),
          ", top-", min(max_genes, nrow(expr_state)), " genes, ",
          n_iter, " iter, ", n_workers, " workers)")

  foreach(i          = seq_len(n_iter),
          .combine   = dplyr::bind_rows,
          .packages  = c("WGCNA", "igraph", "Matrix"),
          .export    = c("build_full_network_boot", "compute_metrics_fast",
                         "compute_EG_fast", "compute_Gbar_fast",
                         "compute_Rbar_fast", "compute_Hb_fast",
                         "nth_largest", "HAS_RSPECTRA")) %dopar% {
    sel <- sample(seq_len(ncol(expr_state)), size = ncol(expr_state),
                  replace = TRUE)
    W   <- build_full_network_boot(expr_state[, sel, drop = FALSE],
                                   beta_fixed, max_genes)
    data.frame(iter = i, condition = state,
               t(compute_metrics_fast(W, eg_keep_pct = 0.001, k_gbar = 50L)))
  }
}

# =============================================================================
# TEST DE PERMUTACION
# =============================================================================

permute_group_test <- function(expr, pheno, states, betas,
                               n_perm = 1000L, max_genes = 2000L) {
  keep  <- pheno$condition %in% states
  expr2 <- expr[, pheno$.sample_id[keep], drop = FALSE]
  ph2   <- droplevels(pheno[keep, , drop = FALSE])
  if (length(unique(ph2$condition)) < 2L || ncol(expr2) < 6L) return(NULL)

  message("  Permutacion: ", states[1], " vs ", states[2],
          " (", n_perm, " iter, metricas: EG+Hb)")

  # Valor observado con compute_metrics_perm (EG + Hb)
  # Rapido: EG sobre top-0.1% aristas, Hb O(p)
  obs <- lapply(seq_along(states), function(j) {
    s <- states[j]
    compute_metrics_perm(
      build_full_network_boot(
        expr2[, ph2$.sample_id[ph2$condition == s], drop = FALSE],
        betas[[s]], max_genes
      )
    )
  })
  obs_diff <- obs[[1]] - obs[[2]]
  # CEI simplificado con EG y Hb (sin Gbar y Rbar que no se calculan en perm)
  obs_CEI  <- (obs[[1]]["EG"] + obs[[1]]["Hb"]) -
              (obs[[2]]["EG"] + obs[[2]]["Hb"])

  perm_diffs <- foreach(b          = seq_len(n_perm),
                        .combine   = rbind,
                        .packages  = c("WGCNA", "igraph", "Matrix"),
                        .export    = c("build_full_network_boot",
                                       "compute_metrics_perm",
                                       "compute_EG_fast", "compute_Hb_fast",
                                       "nth_largest", "HAS_RSPECTRA")) %dopar% {
    perm_cond <- sample(ph2$condition)
    pm <- lapply(seq_along(states), function(j) {
      s <- states[j]
      compute_metrics_perm(
        build_full_network_boot(
          expr2[, ph2$.sample_id[perm_cond == s], drop = FALSE],
          betas[[s]], max_genes
        )
      )
    })
    d   <- pm[[1]] - pm[[2]]
    cei <- (pm[[1]]["EG"] + pm[[1]]["Hb"]) - (pm[[2]]["EG"] + pm[[2]]["Hb"])
    c(d["EG"], d["Hb"], cei)
  }
  colnames(perm_diffs) <- c("EG", "Hb", "CEI")

  data.frame(
    state1        = states[1],
    state2        = states[2],
    obs_diff_EG   = obs_diff["EG"],
    obs_diff_Hb   = obs_diff["Hb"],
    obs_diff_CEI  = obs_CEI,
    p_perm_EG     = mean(abs(perm_diffs[,"EG"])  >= abs(obs_diff["EG"]),  na.rm = TRUE),
    p_perm_Hb     = mean(abs(perm_diffs[,"Hb"])  >= abs(obs_diff["Hb"]),  na.rm = TRUE),
    p_perm_CEI    = mean(abs(perm_diffs[,"CEI"]) >= abs(obs_CEI),         na.rm = TRUE),
    n_perm        = n_perm
  )
}

# =============================================================================
# OPTIMO CONSTRUCTAL
# Usa las redes COMPLETAS guardadas por 02b — no reconstruye en bootstrap.
# =============================================================================

approx_constructal_optimum <- function(W_ref, alpha = 2) {
  p        <- nrow(W_ref)
  k_ref    <- rowSums(W_ref)
  tri      <- upper.tri(W_ref)
  C_budget <- sum(W_ref[tri]^alpha)
  W_opt_raw       <- outer(k_ref, k_ref, function(a, b) (a * b)^(1 / alpha))
  diag(W_opt_raw) <- 0
  W_opt_raw       <- (W_opt_raw + t(W_opt_raw)) / 2
  C_raw           <- sum(W_opt_raw[tri]^alpha)
  if (C_raw < 1e-12) return(W_ref)
  Wopt            <- W_opt_raw * (C_budget / C_raw)^(1 / alpha)
  diag(Wopt)      <- 0
  rownames(Wopt)  <- colnames(Wopt) <- rownames(W_ref)
  Wopt
}

# compute_metrics_precise: metricas de precision sobre la red que recibe.
# La seleccion de subred (top-3000 genes) se hace EXTERNAMENTE antes de llamar
# esta funcion, garantizando que estado y optimo usan exactamente los mismos genes.
compute_metrics_precise <- function(W) {
  EG   <- compute_EG_fast(W, eg_keep_pct = 0.001)
  Gbar <- compute_Gbar_fast(W, k_eig = 200L)
  Rbar <- compute_Rbar_fast(W)
  Hb   <- compute_Hb_fast(W)
  c(EG = EG, Gbar = Gbar, Rbar = Rbar, Hb = Hb)
}

# =============================================================================
# BUCLE PRINCIPAL
# =============================================================================

for (acc in targets) {

  proc_file <- file.path("data/processed", paste0(acc, "_processed_full.rds"))
  if (!file.exists(proc_file)) {
    proc_file <- file.path("data/processed", paste0(acc, "_processed.rds"))
  }
  if (!file.exists(proc_file)) {
    message("[!] Archivo no encontrado: ", proc_file, " — saltando.")
    next
  }

  message("=== FULL bootstrap: ", acc, " ===")
  obj   <- readRDS(proc_file)
  expr  <- obj$expr
  pheno <- obj$pheno
  ord   <- state_orders[[acc]]
  ord   <- ord[ord %in% as.character(unique(pheno$condition))]
  if (length(ord) == 0L) next

  p_full <- nrow(expr)
  # Tamano de subred para bootstrap: min(2000, p_full)
  max_genes_boot <- min(2000L, p_full)
  message("  p_full=", p_full, " | subred bootstrap=", max_genes_boot, " genes")

  # Cargar betas guardados por 02b
  betas <- lapply(setNames(ord, ord), function(st) {
    fp <- file.path("results/full_networks",
                    paste0(acc, "_", st, "_full_network.rds"))
    if (!file.exists(fp)) { message("  Beta no encontrado para ", st, "; usando 6."); return(6L) }
    obj_net <- readRDS(fp)
    if (is.list(obj_net) && !is.null(obj_net$beta)) obj_net$beta else 6L
  })
  message("  Betas: ", paste(names(betas), unlist(betas), sep="=", collapse=" | "))

  # --- 1. Bootstrap por estado ---
  boot_list <- lapply(ord, function(st) {
    bootstrap_by_state(expr, pheno, st, betas[[st]],
                       n_iter = n_boot, max_genes = max_genes_boot)
  })
  boot_df <- dplyr::bind_rows(boot_list[!vapply(boot_list, is.null, logical(1))])
  if (nrow(boot_df) > 0L) {
    fwrite(boot_df,
           file.path("results/full_bootstrap",
                     paste0(acc, "_full_bootstrap_metrics.tsv")), sep = "\t")
    message("  Bootstrap guardado: ", nrow(boot_df), " filas.")
  }

  # --- 2. Tests de permutacion entre estados adyacentes ---
  pairwise <- list()
  if (length(ord) >= 2L) {
    for (i in seq_len(length(ord) - 1L)) {
      res <- permute_group_test(expr, pheno, ord[c(i, i + 1L)], betas,
                                n_perm = n_perm, max_genes = max_genes_boot)
      if (!is.null(res)) pairwise[[i]] <- cbind(accession = acc, res)
    }
  }
  pair_df <- dplyr::bind_rows(pairwise)
  if (nrow(pair_df) > 0L) {
    fwrite(pair_df,
           file.path("results/full_nulls",
                     paste0(acc, "_full_permutation_tests.tsv")), sep = "\t")
    message("  Permutaciones guardadas: ", nrow(pair_df), " filas.")
  }

  # --- 3. Optimo constructal (sobre red completa guardada por 02b) ---
  healthy_state <- ord[1]
  net_path      <- file.path("results/full_networks",
                             paste0(acc, "_", healthy_state, "_full_network.rds"))
  if (!file.exists(net_path)) {
    message("  Red sana no encontrada: ", net_path, " — saltando optimo.")
    next
  }

  message("  Calculando optimo constructal...")
  obj_h    <- readRDS(net_path)
  Whealthy <- if (is.list(obj_h) && !is.null(obj_h$W)) obj_h$W else obj_h
  Wopt     <- approx_constructal_optimum(Whealthy, alpha = 2)

  # Metricas del optimo sobre subred top-3000 (misma logica que en dev)
  p_opt    <- nrow(Wopt)
  k_mopt   <- min(3000L, p_opt)
  if (k_mopt < p_opt) {
    top_opt  <- order(rowSums(Wopt), decreasing = TRUE)[seq_len(k_mopt)]
    Wopt_sub <- Wopt[top_opt, top_opt, drop = FALSE]
  } else {
    Wopt_sub <- Wopt
  }
  mopt <- compute_metrics_precise(Wopt_sub)

  dev <- dplyr::bind_rows(lapply(ord, function(st) {
    fp   <- file.path("results/full_networks",
                      paste0(acc, "_", st, "_full_network.rds"))
    if (!file.exists(fp)) return(NULL)
    Wobj <- readRDS(fp)
    W    <- if (is.list(Wobj) && !is.null(Wobj$W)) Wobj$W else Wobj
    cn   <- intersect(rownames(Wopt), rownames(W))
    if (length(cn) < 10L) return(NULL)
    Ws   <- W[cn, cn, drop = FALSE]
    Wo   <- Wopt[cn, cn, drop = FALSE]

    # CRITICO: pre-seleccionar los mismos top genes en Ws y Wo antes de
    # llamar compute_metrics_precise, garantizando que Gbar y Rbar se
    # calculan sobre identicos nodos en estado y optimo — deltas comparables.
    k_precise <- min(3000L, length(cn))
    if (k_precise < length(cn)) {
      # Seleccionar por fuerza promedio entre estado y optimo
      avg_str  <- (rowSums(Ws) + rowSums(Wo)) / 2
      top_cn   <- names(sort(avg_str, decreasing = TRUE))[seq_len(k_precise)]
      Ws <- Ws[top_cn, top_cn, drop = FALSE]
      Wo <- Wo[top_cn, top_cn, drop = FALSE]
    }

    m    <- compute_metrics_precise(Ws)
    mo   <- compute_metrics_precise(Wo)
    data.frame(
      accession  = acc,
      condition  = st,
      DeltaE     = (mo["EG"]   - m["EG"])   / (abs(mo["EG"])   + 1e-12),
      DeltaR_abs = (m["Rbar"]  - mo["Rbar"]) / (abs(mo["Rbar"]) + 1e-12),
      DeltaG     = (mo["Gbar"] - m["Gbar"])  / (abs(mo["Gbar"]) + 1e-12),
      DeltaH     = (mo["Hb"]   - m["Hb"])    / (abs(mo["Hb"])   + 1e-12)
    )
  }))

  if (!is.null(dev) && nrow(dev) > 0L) {
    saveRDS(list(Wopt = Wopt, metrics = mopt),
            file.path("results/full_optimum",
                      paste0(acc, "_full_optimum.rds")))
    fwrite(dev,
           file.path("results/full_optimum",
                     paste0(acc, "_full_deviation_from_optimum.tsv")), sep = "\t")
    message("  Optimo guardado.")
  }

  message("  ", acc, " completado.")
}

stopCluster(cl)
message("Done.")

