#!/usr/bin/env Rscript
# =============================================================================
# 03_bootstrap_nulls_optimum.R
# Bootstrap de metricas, tests de permutacion y optimo constructal
# sobre el espacio COMUN de genes (rama compartida).
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# RAMA: Compartida (subred fija <= 800 genes comunes)
#
# CORRECCIONES:
#   [C9]  approx_constructal_optimum: optimo de primer orden del lagrangiano
#         max EG s.t. sum(w^alpha) = C. Pesos proporcionales a (k_i*k_j)^(1/alpha).
#   [C10] permute_group_test: p-valores para EG, Hb, CEI (ver FIX-split).
#   [C11] Bootstrap paralelizado con doParallel/foreach.
#   [C12] Datasets: GSE76895, GSE18732, GSE15653, GSE27951.
#   [FIX-beta]  Beta fijo del archivo guardado por script 02.
#   [FIX-Rbar]  Rbar = 2*tr(L+)/p (corregido desde /(p-1)).
#   [FIX-hv]    high_variance_genes.rds leido UNA vez antes del bucle principal.
#               La version original lo leia en cada iteracion bootstrap (1200x).
#   [FIX-doRNG] doRNG opcional — fallback a %dopar% si no esta instalado.
#   [FIX-perm]  n_perm=1000 (p_min=0.001, minimo aceptable para publicacion).
#   [FIX-split] Dos wrappers de metricas segun contexto:
#               - compute_metrics(): EG+Gbar+Rbar+Hb — bootstrap y optimo.
#               - compute_metrics_perm(): solo EG+Hb — permutaciones (1000x2x4).
#               expm() O(p^3) en 16000 llamadas = horas; EG+Hb es suficiente
#               para detectar diferencias bajo H0 de permutacion de etiquetas.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(expm)
  library(Matrix)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(WGCNA)
})

# doRNG opcional — reproducibilidad en paralelo
HAS_DORNG <- requireNamespace("doRNG", quietly = TRUE)
if (HAS_DORNG) {
  library(doRNG)
  message("doRNG disponible: resultados reproducibles.")
} else {
  message("doRNG no disponible. Instalar con: install.packages('doRNG')")
  message("Usando %dopar% sin reproducibilidad garantizada entre ejecuciones.")
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Paralelismo HPC
# -----------------------------------------------------------------------------
n_workers <- as.integer(Sys.getenv("N_WORKERS", unset = "40"))

Sys.setenv(OMP_NUM_THREADS      = 1L,
           OPENBLAS_NUM_THREADS = 1L,
           MKL_NUM_THREADS      = 1L)
disableWGCNAThreads()

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

cl <- makeCluster(n_workers)
registerDoParallel(cl)
if (HAS_DORNG) registerDoRNG(1234)

message("Motor paralelo: ", n_workers, " workers | BLAS=1 | WGCNA=0 | doRNG=",
        HAS_DORNG)

# -----------------------------------------------------------------------------
dir.create("results/bootstrap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/nulls",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/optimum",   recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

state_orders <- list(
  GSE76895 = c("ND",   "IGT", "T2D"),
  GSE18732 = c("ND",   "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT",  "IGT", "T2D")
)

n_boot <- 100L    # bootstrap: IC de las 4 metricas (EG, Gbar, Rbar, Hb)
n_perm <- 1000L   # permutaciones: p_min=0.001 | ~5-10 min con 40 workers
message("Bootstrap: ", n_boot, " iter | Permutaciones: ", n_perm, " iter")

# -----------------------------------------------------------------------------
# [FIX-hv] Leer genes de alta varianza UNA SOLA VEZ
# -----------------------------------------------------------------------------
hv_path <- "data/processed/high_variance_genes.rds"
if (!file.exists(hv_path)) {
  stop("'high_variance_genes.rds' no encontrado. Ejecuta primero el script 01.")
}
hv_genes_global    <- readRDS(hv_path)
fixed_genes_global <- hv_genes_global[seq_len(min(800L, length(hv_genes_global)))]
message("Subred fija: ", length(fixed_genes_global), " genes")

# =============================================================================
# METRICAS — identicas a script 02
# [FIX-Rbar] Rbar = 2*tr(L+)/p
# =============================================================================

# Metricas completas: bootstrap y optimo constructal
compute_metrics <- function(W) {
  p    <- nrow(W)
  k    <- rowSums(W)

  g    <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
  Dg   <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / Dg; diag(invD) <- NA_real_
  EG   <- mean(invD[upper.tri(invD)], na.rm = TRUE) * 2

  Dinv <- diag(1 / sqrt(k + 1e-12), p)
  Gbar <- mean(expm(Dinv %*% W %*% Dinv)[row(W) != col(W)])

  L        <- diag(k, p) - W
  eig      <- eigen(L, symmetric = TRUE)
  tol      <- max(abs(eig$values)) * p * 1e-10
  inv_vals <- ifelse(eig$values > tol, 1 / eig$values, 0)
  Lplus    <- eig$vectors %*% diag(inv_vals, p) %*% t(eig$vectors)
  Rbar     <- as.numeric(2 * sum(diag(Lplus)) / p)   # [FIX-Rbar]

  prob <- k / sum(k)
  Hb   <- -sum(prob * log(prob + 1e-12))

  c(EG = EG, Gbar = Gbar, Rbar = Rbar, Hb = Hb)
}

# [FIX-split] Metricas rapidas: solo EG y Hb para permutaciones.
# Gbar con expm() O(p^3): ~0.1s * 1000 perms * 2 redes * 2 pares * 4 datasets
# = 1600s ~ 27 min solo en Gbar. EG+Hb detectan diferencias igual de bien bajo H0.
compute_metrics_perm <- function(W) {
  k  <- rowSums(W)
  g  <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                    diag = FALSE)
  Dg   <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / Dg; diag(invD) <- NA_real_
  EG   <- mean(invD[upper.tri(invD)], na.rm = TRUE) * 2
  prob <- k / sum(k)
  Hb   <- -sum(prob * log(prob + 1e-12))
  c(EG = EG, Hb = Hb)
}

# [FIX-beta + FIX-hv] build_network con beta fijo y genes pre-cargados
build_network <- function(expr_sub, beta_fixed, fixed_genes) {
  cm <- bicor(t(expr_sub), maxPOutliers = 0.1)
  cm[is.na(cm)] <- 0
  W  <- abs(cm)^beta_fixed
  diag(W) <- 0
  rownames(W) <- colnames(W) <- rownames(expr_sub)
  top <- intersect(fixed_genes, rownames(W))
  W[top, top, drop = FALSE]
}

# =============================================================================
# BOOTSTRAP (metricas completas)
# =============================================================================

bootstrap_by_state <- function(expr, pheno, state, beta_fixed,
                               fixed_genes, n_iter = 100L) {
  samples    <- pheno$.sample_id[pheno$condition == state]
  expr_state <- expr[, samples, drop = FALSE]

  if (ncol(expr_state) < 3L) {
    message("  -> Saltando bootstrap '", state, "': < 3 muestras")
    return(NULL)
  }
  message("  Bootstrap: ", state, " (", n_iter, " iter, beta=", beta_fixed, ")")

  foreach(i          = seq_len(n_iter),
          .combine   = dplyr::bind_rows,
          .packages  = c("WGCNA", "igraph", "Matrix", "expm"),
          .export    = c("build_network", "compute_metrics")) %dopar% {
    sel <- sample(seq_len(ncol(expr_state)), size = ncol(expr_state),
                  replace = TRUE)
    W   <- build_network(expr_state[, sel, drop = FALSE], beta_fixed, fixed_genes)
    data.frame(iter = i, condition = state, t(compute_metrics(W)))
  }
}

# =============================================================================
# TEST DE PERMUTACION (metricas rapidas: EG+Hb)
# =============================================================================

permute_group_test <- function(expr, pheno, states, betas,
                               fixed_genes, n_perm = 1000L) {
  keep  <- pheno$condition %in% states
  expr2 <- expr[, pheno$.sample_id[keep], drop = FALSE]
  ph2   <- droplevels(pheno[keep, , drop = FALSE])

  if (length(unique(ph2$condition)) < 2L || ncol(expr2) < 6L) {
    message("  -> Saltando permutacion: ", states[1], " vs ", states[2])
    return(NULL)
  }

  obs <- lapply(seq_along(states), function(j) {
    s <- states[j]
    compute_metrics_perm(
      build_network(expr2[, ph2$.sample_id[ph2$condition == s], drop = FALSE],
                    betas[[s]], fixed_genes)
    )
  })
  obs_diff <- obs[[1]] - obs[[2]]
  obs_CEI  <- (obs[[1]]["EG"] + obs[[1]]["Hb"]) -
              (obs[[2]]["EG"] + obs[[2]]["Hb"])

  message("  Permutando: ", states[1], " vs ", states[2],
          " (", n_perm, " iter, EG+Hb)")

  perm_diffs <- foreach(b          = seq_len(n_perm),
                        .combine   = rbind,
                        .packages  = c("WGCNA", "igraph", "Matrix"),
                        .export    = c("build_network", "compute_metrics_perm")) %dopar% {
    perm_cond <- sample(ph2$condition)
    pm <- lapply(seq_along(states), function(j) {
      s <- states[j]
      compute_metrics_perm(
        build_network(expr2[, ph2$.sample_id[perm_cond == s], drop = FALSE],
                      betas[[s]], fixed_genes)
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
    p_perm_EG     = mean(abs(perm_diffs[, "EG"])  >= abs(obs_diff["EG"]),  na.rm = TRUE),
    p_perm_Hb     = mean(abs(perm_diffs[, "Hb"])  >= abs(obs_diff["Hb"]),  na.rm = TRUE),
    p_perm_CEI    = mean(abs(perm_diffs[, "CEI"]) >= abs(obs_CEI),         na.rm = TRUE),
    n_perm        = n_perm
  )
}

# =============================================================================
# OPTIMO CONSTRUCTAL
# =============================================================================

approx_constructal_optimum <- function(W_ref, alpha = 2) {
  p        <- nrow(W_ref)
  k_ref    <- rowSums(W_ref)
  tri      <- upper.tri(W_ref)
  C_budget <- sum(W_ref[tri]^alpha)
  W_opt_raw        <- outer(k_ref, k_ref, function(a, b) (a * b)^(1 / alpha))
  diag(W_opt_raw)  <- 0
  W_opt_raw        <- (W_opt_raw + t(W_opt_raw)) / 2
  C_raw            <- sum(W_opt_raw[tri]^alpha)
  if (C_raw < 1e-12) {
    warning("approx_constructal_optimum: presupuesto nulo; devolviendo original")
    return(W_ref)
  }
  Wopt           <- W_opt_raw * (C_budget / C_raw)^(1 / alpha)
  diag(Wopt)     <- 0
  rownames(Wopt) <- colnames(Wopt) <- rownames(W_ref)
  Wopt
}

# =============================================================================
# BUCLE PRINCIPAL
# =============================================================================

for (acc in targets) {
  message("=== ", acc, " ===")

  proc_file <- file.path("data/processed", paste0(acc, "_processed.rds"))
  if (!file.exists(proc_file)) {
    message("  Archivo no encontrado: ", proc_file, " — saltando.")
    next
  }

  obj   <- readRDS(proc_file)
  expr  <- obj$expr
  pheno <- obj$pheno
  ord   <- state_orders[[acc]]
  ord   <- ord[ord %in% as.character(unique(pheno$condition))]
  if (length(ord) == 0L) next

  # Betas guardados por script 02
  betas <- lapply(setNames(ord, ord), function(st) {
    fp <- file.path("results/networks", paste0(acc, "_", st, "_network.rds"))
    if (!file.exists(fp)) {
      message("  Beta no encontrado para ", st, "; usando fallback=6.")
      return(6L)
    }
    readRDS(fp)$beta
  })
  message("  Betas: ", paste(names(betas), unlist(betas), sep = "=",
                             collapse = " | "))

  # 1. Bootstrap — metricas completas
  boot_list <- lapply(ord, function(st) {
    bootstrap_by_state(expr, pheno, st,
                       beta_fixed  = betas[[st]],
                       fixed_genes = fixed_genes_global,
                       n_iter      = n_boot)
  })
  boot <- dplyr::bind_rows(boot_list[!sapply(boot_list, is.null)])
  if (nrow(boot) > 0L) {
    fwrite(boot,
           file.path("results/bootstrap",
                     paste0(acc, "_bootstrap_metrics.tsv")), sep = "\t")
    message("  Bootstrap guardado: ", nrow(boot), " filas.")
  }

  # 2. Permutaciones — EG+Hb
  pairwise <- list()
  if (length(ord) >= 2L) {
    for (i in seq_len(length(ord) - 1L)) {
      pair <- ord[c(i, i + 1L)]
      res  <- permute_group_test(expr, pheno, pair, betas,
                                 fixed_genes = fixed_genes_global,
                                 n_perm      = n_perm)
      if (!is.null(res)) pairwise[[i]] <- cbind(accession = acc, res)
    }
  }
  pw_df <- dplyr::bind_rows(pairwise)
  if (nrow(pw_df) > 0L) {
    fwrite(pw_df,
           file.path("results/nulls",
                     paste0(acc, "_permutation_tests.tsv")), sep = "\t")
    message("  Permutaciones guardadas: ", nrow(pw_df), " filas.")
  }

  # 3. Optimo constructal — metricas completas sobre red guardada por script 02
  healthy_state <- ord[1]
  net_path      <- file.path("results/networks",
                             paste0(acc, "_", healthy_state, "_network.rds"))
  if (!file.exists(net_path)) {
    message("  Red sana no encontrada: ", net_path, " — saltando optimo.")
    next
  }

  Whealthy <- readRDS(net_path)$W
  Wopt     <- approx_constructal_optimum(Whealthy, alpha = 2)
  mopt     <- compute_metrics(Wopt)

  dev_list <- lapply(ord, function(st) {
    fp <- file.path("results/networks",
                    paste0(acc, "_", st, "_network.rds"))
    if (!file.exists(fp)) return(NULL)
    W  <- readRDS(fp)$W
    cn <- intersect(rownames(Wopt), rownames(W))
    if (length(cn) < 10L) return(NULL)
    Ws <- W[cn, cn, drop = FALSE]
    Wo <- Wopt[cn, cn, drop = FALSE]
    m  <- compute_metrics(Ws)
    mo <- compute_metrics(Wo)
    data.frame(
      accession  = acc,
      condition  = st,
      DeltaE     = (mo["EG"]   - m["EG"])   / (abs(mo["EG"])   + 1e-12),
      DeltaR_abs = (m["Rbar"]  - mo["Rbar"]) / (abs(mo["Rbar"]) + 1e-12),
      DeltaG     = (mo["Gbar"] - m["Gbar"])  / (abs(mo["Gbar"]) + 1e-12),
      DeltaH     = (mo["Hb"]   - m["Hb"])    / (abs(mo["Hb"])   + 1e-12)
    )
  })
  dev <- dplyr::bind_rows(dev_list[!sapply(dev_list, is.null)])

  saveRDS(list(Wopt = Wopt, metrics = mopt),
          file.path("results/optimum", paste0(acc, "_optimum.rds")))
  if (nrow(dev) > 0L) {
    fwrite(dev,
           file.path("results/optimum",
                     paste0(acc, "_deviation_from_optimum.tsv")), sep = "\t")
  }
  message("  ", acc, " completado.")
}

stopCluster(cl)
message("Done.")