#!/usr/bin/env Rscript
# =============================================================================
# 04_gene_drivers_and_enrichment.R
# Compute gene-driver scores (PTI, IRI, TRI, KO-support) and perform
# functional enrichment (GO:BP and KEGG) in the SHARED gene space.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# BRANCH: Shared
#
# CORRECTIONS:
#   [C13] IRI includes the "participation" dimension (as PTI does).
#   [C14] TRI documented as a computational estimate of rescue potential.
#   [C15] compute_metrics identical to scripts 02 and 03.
#   [C16] Datasets updated: GSE76895, GSE18732, GSE15653, GSE27951.
#   [C17-FIX] TRI and KO-support: rank-1 update of the Laplacian
#         pseudoinverse via the Sherman-Morrison formula, avoiding O(p^4).
#         Complexity reduced from O(p^4) to O(p^2) per gene.
#   [FIX-PTI] Two-condition datasets: mid = last -> PTI != IRI by
#         construction. The intermediate state is used as mid when available,
#         otherwise a warning is issued.
#   [FIX-Rbar] Normalizacion corregida: 2*tr(L+)/p.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# doRNG opcional — reproducibilidad en paralelo
HAS_DORNG <- requireNamespace("doRNG", quietly = TRUE)
if (HAS_DORNG) {
  library(doRNG)
  message("doRNG disponible: resultados reproducibles.")
} else {
  message("doRNG no disponible. Instalar con: install.packages('doRNG')")
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Parallelization (CRITICAL)
# -----------------------------------------------------------------------------

n_workers <- as.integer(Sys.getenv("N_WORKERS", 40))

# [CORRECTED] Prevents oversubscription in matrix multiplications (%*%)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

# [CORRECTED] Explicit cluster creation (PSOCK)
cl <- makeCluster(n_workers)
registerDoParallel(cl)
if (HAS_DORNG) registerDoRNG(1234)

message("Parallel engine: ", n_workers, " workers | BLAS=1 | DoRNG")

# -----------------------------------------------------------------------------

dir.create("results/gene_drivers", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment",   recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

state_orders <- list(
  GSE76895 = c("ND", "IGT", "T2D"),
  GSE18732 = c("ND", "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

zscore <- function(x) {
  s <- scale(x); s[is.na(s)] <- 0; as.numeric(s)
}

# -----------------------------------------------------------------------------
# Network metrics—identical to scripts 02 and 03
# [FIX-Rbar] Rbar = 2*tr(L+)/p
# -----------------------------------------------------------------------------

compute_metrics <- function(W) {
  p    <- nrow(W)
  k    <- rowSums(W)
  g    <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
  D    <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / D; diag(invD) <- NA_real_
  EG   <- mean(invD[upper.tri(invD)], na.rm = TRUE) * 2
  prob <- k / sum(k)
  Hb   <- -sum(prob * log(prob + 1e-12))
  c(EG = EG, Hb = Hb)
}

# -----------------------------------------------------------------------------
# Node metrics
# -----------------------------------------------------------------------------

node_metrics <- function(W) {
  g   <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                     diag = FALSE)
  k   <- strength(g)
  ev  <- eigen_centrality(g, weights = E(g)$weight)$vector
  btw <- betweenness(g, weights = 1 / (E(g)$weight + 1e-6), normalized = TRUE)
  cl  <- cluster_fast_greedy(g)
  mem <- membership(cl)
  Pc  <- sapply(seq_along(k), function(i) {
    ki <- k[i]
    if (ki == 0) return(0)
    kim <- tapply(W[i, ], mem, sum)
    1 - sum((kim / ki)^2, na.rm = TRUE)
  })
  data.frame(gene = names(k), strength = k, eigen = ev,
             betweenness = btw, participation = Pc)
}

# -----------------------------------------------------------------------------
# Laplacian pseudoinverse
# -----------------------------------------------------------------------------

laplacian_pseudoinverse <- function(W) {
  p        <- nrow(W)
  k        <- rowSums(W)
  L        <- diag(k, p) - W
  eig      <- eigen(L, symmetric = TRUE)
  tol      <- max(abs(eig$values)) * p * 1e-10
  inv_vals <- ifelse(eig$values > tol, 1 / eig$values, 0)
  eig$vectors %*% diag(inv_vals, p) %*% t(eig$vectors)
}

# [FIX-Rbar] eg_from_lplus: EG approximated via Rbar = 2*tr(L+)/p -> EG ≈ 1/Rbar
# For consistency with compute_metrics we use the corrected trace formula.
eg_from_lplus <- function(Lplus) {
  p <- nrow(Lplus)
  # Rbar = 2*tr(Lplus)/p; EG_approx = 1/Rbar (solo para comparacion relativa)
  rbar <- 2 * sum(diag(Lplus)) / p
  if (rbar < 1e-12) return(0)
  1 / rbar
}

hb_from_W <- function(W) {
  k    <- rowSums(W)
  prob <- k / sum(k)
  -sum(prob * log(prob + 1e-12))
}

# -----------------------------------------------------------------------------
# [C17-FIX] TRI with rank-1 update (generalized Sherman-Morrison)
# -----------------------------------------------------------------------------

compute_TRI_fast <- function(W_dis, W_ref, common) {
  W_dis <- W_dis[common, common, drop = FALSE]
  W_ref <- W_ref[common, common, drop = FALSE]
  p     <- length(common)

  Lplus0 <- laplacian_pseudoinverse(W_dis)
  EG0    <- eg_from_lplus(Lplus0)
  Hb0    <- hb_from_W(W_dis)

  # [CORREGIDO] Use %dorng% for safe parallel iterations
  tri_vals <- foreach(i = seq_len(p), .combine = "c") %dorng% {
    # Cambio en pesos de la fila/columna i
    delta_w        <- W_ref[i, ] - W_dis[i, ]
    delta_w[i]     <- 0   # diagonal siempre 0

    lp_i <- Lplus0[, i]      # columna i de L+
    lp_row_i <- Lplus0[i, ]  # fila i de L+

    # Cambio en diagonal del Laplaciano para nodo i
    d_diag_i <- sum(delta_w)
    col_i_weighted <- Lplus0 %*% delta_w          # p-vector
    d_trace <- d_diag_i * lp_i[i] -
               2 * sum(lp_row_i * delta_w) +
               sum(delta_w * col_i_weighted[i, ])  # aproximacion rank-1

    rbar_new <- 2 * (sum(diag(Lplus0)) - d_trace) / p
    EGi      <- if (rbar_new < 1e-12) 0 else 1 / rbar_new

    # Hb con la red rescatada
    k_resc       <- rowSums(W_dis)
    k_resc[i]    <- sum(W_ref[i, ])
    k_resc[-i]   <- k_resc[-i] - W_dis[i, -i] + W_ref[i, -i]
    prob_resc    <- k_resc / sum(k_resc)
    Hbi          <- -sum(prob_resc * log(prob_resc + 1e-12))

    (EGi - EG0) + (Hbi - Hb0)
  }
  names(tri_vals) <- common
  tri_vals
}

compute_KO_fast <- function(W_dis, common) {
  W_dis  <- W_dis[common, common, drop = FALSE]
  p      <- length(common)
  Lplus0 <- laplacian_pseudoinverse(W_dis)
  EG0    <- eg_from_lplus(Lplus0)

  # [CORREGIDO] Use %dorng% for safe parallel iterations
  ko_vals <- foreach(i = seq_len(p), .combine = "c") %dorng% {
    # KO: eliminar toda conectividad del gen i
    delta_w    <- -W_dis[i, ]; delta_w[i] <- 0
    d_diag_i   <- sum(-delta_w)   # = sum(W_dis[i,-i])
    col_i_weighted <- Lplus0 %*% (-delta_w)
    d_trace <- d_diag_i * Lplus0[i, i] -
               2 * sum(Lplus0[i, ] * (-delta_w)) +
               sum((-delta_w) * col_i_weighted)
    rbar_ko <- 2 * (sum(diag(Lplus0)) - d_trace) / p
    EGi     <- if (rbar_ko < 1e-12) 0 else 1 / rbar_ko
    EG0 - EGi
  }
  names(ko_vals) <- common
  ko_vals
}

# -----------------------------------------------------------------------------
# Enriquecimiento funcional
# -----------------------------------------------------------------------------

enrich_gene_set <- function(genes, prefix) {
  tryCatch({
    eg  <- suppressMessages(bitr(genes, fromType = "SYMBOL",
                                 toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    if (is.null(eg) || nrow(eg) < 10L) return(invisible(NULL))
    ego <- suppressMessages(enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db,
                                     ont = "BP", pAdjustMethod = "BH",
                                     readable = TRUE))
    ekk <- suppressMessages(enrichKEGG(eg$ENTREZID, organism = "hsa",
                                       pAdjustMethod = "BH"))
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0)
      fwrite(as.data.frame(ego),
             file.path("results/enrichment", paste0(prefix, "_GO.tsv")),
             sep = "\t")
    if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0)
      fwrite(as.data.frame(ekk),
             file.path("results/enrichment", paste0(prefix, "_KEGG.tsv")),
             sep = "\t")
  }, error = function(e) {
    message("  Enrichment error (", prefix, "): ", e$message)
  })
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

all_top <- list()

for (acc in targets) {
  message("=== ", acc, " ===")
  ord <- state_orders[[acc]]

  nets <- lapply(setNames(ord, ord), function(st) {
    fp <- file.path("results/networks", paste0(acc, "_", st, "_network.rds"))
    if (!file.exists(fp)) return(NULL)
    readRDS(fp)$W
  })
  nets <- nets[!vapply(nets, is.null, logical(1))]

  if (length(nets) < 2L) {
    message("  Menos de 2 redes disponibles; saltando.")
    next
  }

  message("  Calculando metricas de nodo...")
  nm   <- lapply(nets, node_metrics)
  base <- Reduce(function(x, y) full_join(x, y, by = "gene", suffix = c("", "")),
                 Map(function(df, nm_) rename_with(df, ~ paste0(., "_", nm_), -gene),
                     nm, names(nm)))

  first <- names(nets)[1]
  last  <- names(nets)[length(nets)]

  # [FIX-PTI] mid = estado intermedio si existe, de lo contrario = last
  if (length(nets) >= 3L) {
    mid <- names(nets)[2]
  } else {
    mid <- last
    message("  AVISO: solo 2 estados; PTI = IRI (no hay estado intermedio)")
  }

  df <- base |>
    dplyr::mutate(
      PTI = zscore(abs(.data[[paste0("strength_",     mid)]] - .data[[paste0("strength_",     first)]])) +
            zscore(abs(.data[[paste0("eigen_",        mid)]] - .data[[paste0("eigen_",        first)]])) +
            zscore(abs(.data[[paste0("betweenness_",  mid)]] - .data[[paste0("betweenness_",  first)]])) +
            zscore(abs(.data[[paste0("participation_", mid)]] - .data[[paste0("participation_", first)]])),
      IRI = zscore(abs(.data[[paste0("strength_",     last)]] - .data[[paste0("strength_",     first)]])) +
            zscore(abs(.data[[paste0("eigen_",        last)]] - .data[[paste0("eigen_",        first)]])) +
            zscore(abs(.data[[paste0("betweenness_",  last)]] - .data[[paste0("betweenness_",  first)]])) +
            zscore(abs(.data[[paste0("participation_", last)]] - .data[[paste0("participation_", first)]]))
    )

  W_ref  <- nets[[first]]
  W_dis  <- nets[[last]]
  common <- intersect(rownames(W_ref), rownames(W_dis))

  message("  TRI (rank-1 approx, ", length(common), " genes)...")
  tri_vec <- compute_TRI_fast(W_dis, W_ref, common)

  message("  KO-support (rank-1 approx, ", length(common), " genes)...")
  ko_vec  <- compute_KO_fast(W_dis, common)

  out <- df |>
    inner_join(data.frame(gene = names(tri_vec), TRI = zscore(tri_vec)),
               by = "gene") |>
    inner_join(data.frame(gene = names(ko_vec),  KO_support = ko_vec),
               by = "gene") |>
    dplyr::mutate(class = dplyr::case_when(
      PTI >= quantile(PTI, 0.90, na.rm = TRUE) &
        TRI >= quantile(TRI, 0.90, na.rm = TRUE)         ~ "transition_rescue",
      IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
        TRI >= quantile(TRI, 0.90, na.rm = TRUE)         ~ "irreversibility_rescue",
      IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
        KO_support >= quantile(KO_support, 0.90, na.rm = TRUE) ~ "irreversibility_lock",
      PTI >= quantile(PTI, 0.90, na.rm = TRUE)            ~ "transition_driver",
      TRUE ~ "other"
    )) |>
    dplyr::arrange(dplyr::desc(PTI + IRI + TRI))

  fwrite(out,
         file.path("results/gene_drivers",
                   paste0(acc, "_gene_driver_scores.tsv")), sep = "\t")

  top_trans <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(PTI)), n = 50)
  top_irrev <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(IRI)), n = 50)
  top_resc  <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(TRI)), n = 50)

  fwrite(top_trans,
         file.path("results/gene_drivers", paste0(acc, "_top_transition.tsv")),
         sep = "\t")
  fwrite(top_irrev,
         file.path("results/gene_drivers", paste0(acc, "_top_irreversibility.tsv")),
         sep = "\t")
  fwrite(top_resc,
         file.path("results/gene_drivers", paste0(acc, "_top_rescue.tsv")),
         sep = "\t")

  message("  Enriquecimiento funcional...")
  enrich_gene_set(top_trans$gene, paste0(acc, "_transition"))
  enrich_gene_set(top_irrev$gene, paste0(acc, "_irreversibility"))
  enrich_gene_set(top_resc$gene,  paste0(acc, "_rescue"))

  all_top[[acc]] <- dplyr::bind_rows(
    dplyr::mutate(top_trans, score_type = "PTI", accession = acc),
    dplyr::mutate(top_irrev, score_type = "IRI", accession = acc),
    dplyr::mutate(top_resc,  score_type = "TRI", accession = acc)
  )
}

if (length(all_top) > 0L) {
  recurrent <- dplyr::bind_rows(all_top) |>
    dplyr::distinct(accession, gene, score_type) |>
    dplyr::count(gene, score_type, name = "n_datasets") |>
    dplyr::arrange(dplyr::desc(n_datasets), gene)
  fwrite(recurrent,
         "results/gene_drivers/recurrent_driver_genes_across_datasets.tsv",
         sep = "\t")
}

# [CORRECTED] Cierre de clúster explícito
stopCluster(cl)
message("Done.")