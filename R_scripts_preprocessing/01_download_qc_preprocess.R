#!/usr/bin/env Rscript
# =============================================================================
# 01_download_qc_preprocess.R
# Descarga datasets GEO, mapea condiciones, colapsa duplicados y filtra genes
# de alta varianza.
#
# APPLIED CORRECTIONS:
#   [C1] collapse_duplicate_genes: removed the native pipe with placeholder '.'
#        incompatible across R versions. Replaced with explicit
#        assignment.
#   [C2] Conditional log2 transformation: checks the maximum range
#        before applying log2(x + 0.01) to avoid double transformation of
#        datasets already returned normalized by GEOquery.
#   [C3] GSE18732: "NGT" relabeled as "ND" to standardize the
#        healthy reference across all datasets.
#   [C4] Global variance filter: the Q50 threshold (median) is computed on the
#        mean variance of genes shared across all datasets, ensuring that
#        nodes are comparable between datasets. The shared-gene intersection
#        (common_genes.rds) is exported for use in script 02.
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(WGCNA)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
})

options(stringsAsFactors = FALSE)
#allowWGCNAThreads()

dir.create("data/raw",       recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc",     recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
# The "metabolic quartet": pancreas, muscle, liver, and adipose tissue
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

safe_coldata <- function(eset) {
  pheno            <- pData(eset)
  pheno$.sample_id <- rownames(pheno)
  pheno
}

# Explicit rowname assignment and variance-based probe collapse
collapse_duplicate_genes <- function(expr, symbols) {
  keep    <- !is.na(symbols) & symbols != "" & !duplicated(rownames(expr))
  expr    <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  idx     <- split(seq_along(symbols), symbols)
  chosen  <- vapply(idx, function(ix) {
    if (length(ix) == 1L) return(ix)
    vars <- apply(expr[ix, , drop = FALSE], 1, var, na.rm = TRUE)
    ix[which.max(vars)]
  }, integer(1))
  out           <- expr[chosen, , drop = FALSE]
  rownames(out) <- symbols[chosen]
  out
}

infer_symbol_column <- function(fdata) {
  cand <- c("Gene Symbol","GENE_SYMBOL","Symbol","Gene symbol",
            "Gene.symbol","SYMBOL","gene_assignment")
  hit  <- intersect(cand, colnames(fdata))
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

extract_gene_symbols <- function(eset, accession) {
  fdata   <- fData(eset)
  sym_col <- infer_symbol_column(fdata)
  
  # 1. Standard attempt: if a standard column is available
  if (!is.na(sym_col)) {
    x <- as.character(fdata[[sym_col]])
    if (sym_col == "gene_assignment") x <- sub(".* // ([^ ]+) //.*", "\\1", x)
    x <- str_replace_all(x, "///.*$", "")
    x <- str_trim(x)
    return(x)
  }
  
  probe_ids <- rownames(exprs(eset))
  
  # 2. Plan B.1: detect Ensembl-based custom CDFs (the GSE18732 case)
  if (any(grepl("ENST|ENSG", head(probe_ids, 100)) | grepl("ENST|ENSG", tail(probe_ids, 100)))) {
    message("  ->   -> Detected Ensembl custom CDF in ", accession, ". mapping with org.Hs.eg.db...")
    
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      warning("Missing package 'org.Hs.eg.db'. Returning original IDs.")
      return(probe_ids)
    }
    
    # Clean Brainarray suffixes (e.g., ENST00000123_at -> ENST00000123)
    clean_ids <- sub("_at$", "", probe_ids)
    clean_ids <- sub("_st$", "", clean_ids)
    
    # Determine whether IDs are transcripts (ENST) or genes (ENSG)
    key_type <- ifelse(any(grepl("ENST", clean_ids)), "ENSEMBLTRANS", "ENSEMBL")
    
    mapped_symbols <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys = clean_ids,
                            column = "SYMBOL",
                            keytype = key_type,
                            multiVals = "first")
    )
    return(unname(mapped_symbols))
  }
  
  # 3. Plan B.2: classic Affymetrix without annotation
  message("  ->   -> No annotation column found in ", accession, ". trying mapping with hthgu133a.db...")
  if (!requireNamespace("hthgu133a.db", quietly = TRUE)) {
    return(probe_ids)
  }
  
  mapped_symbols <- suppressMessages(
    AnnotationDbi::mapIds(hthgu133a.db::hthgu133a.db,
                          keys = probe_ids,
                          column = "SYMBOL",
                          keytype = "PROBEID",
                          multiVals = "first")
  )
  return(unname(mapped_symbols))
}

# Ad hoc mapping for the four selected datasets
map_conditions <- function(accession, pheno) {
  # Texto colapsado general (seguro para GSE76895 y GSE27951)
  txt   <- apply(pheno, 1, function(x) paste(x, collapse = " | "))
  txt_l <- str_to_lower(txt)
  
  # Extraer la columna title en minúsculas de forma segura
  title_l <- if ("title" %in% colnames(pheno)) str_to_lower(pheno$title) else txt_l
  
  condition <- rep(NA_character_, length(txt_l))

  if (accession == "GSE76895") {
    condition[str_detect(txt_l, "\\| nd \\|")]  <- "ND"
    condition[str_detect(txt_l, "\\| igt \\|")] <- "IGT"
    condition[str_detect(txt_l, "\\| t2d \\|")] <- "T2D"
    condition[str_detect(txt_l, "\\| t3cd \\|")] <- "T3cD"

  } else if (accession == "GSE18732") {
    # Usamos SOLO la columna title, que es súper limpia
    condition[str_detect(title_l, "normal")]            <- "ND"
    condition[str_detect(title_l, "glucoseintolerant")] <- "IGT"
    condition[str_detect(title_l, "diabetic")]          <- "T2D"

  } else if (accession == "GSE15653") {
    # Usamos SOLO la columna title
    condition[str_detect(title_l, "lean")]       <- "Lean"
    condition[str_detect(title_l, "obese_nodm")] <- "Obese_noT2D"
    condition[str_detect(title_l, "obese_dm")]   <- "Obese_T2D"
    
  } else if (accession == "GSE27951") {
    condition[str_detect(txt_l, "\\b(ngt|normal)\\b")] <- "NGT"
    condition[str_detect(txt_l, "\\b(igt|impaired)\\b")] <- "IGT"
    condition[str_detect(txt_l, "\\b(type 2|t2d|dm)\\b")] <- "T2D"
  }

  data.frame(sample_id = pheno$.sample_id, condition = condition)
}

# Fail-safe scale heuristic
needs_log2 <- function(expr) {
  if (length(expr) == 0 || all(is.na(expr))) return(FALSE)
  max_val <- suppressWarnings(max(expr, na.rm = TRUE))
  if (is.infinite(max_val)) return(FALSE)
  max_val > 30
}

# -----------------------------------------------------------------------------
# First pass: download and store raw objects
# -----------------------------------------------------------------------------

preprocess_dataset <- function(accession) {
  message("Processing ", accession)
  gse <- getGEO(accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
  
  # Smart matrix selection (protection for GSE27951 and similar cases)
  if (length(gse) > 1) {
    message("  ->   -> Multiple matrices detected in ", accession, ". searching for GPL570...")
    idx <- grep("GPL570", names(gse))
    if (length(idx) > 0) {
      eset <- gse[[idx[1]]]
    } else {
      eset <- gse[[1]]
    }
  } else {
    eset <- gse[[1]]
  }
  
  message("  ->   -> Selected platform: ", annotation(eset))
  
  expr  <- exprs(eset)
  pheno <- safe_coldata(eset)
  cond  <- map_conditions(accession, pheno)
  pheno <- left_join(pheno, cond, by = c(".sample_id" = "sample_id"))

  symbols <- extract_gene_symbols(eset, accession)
  expr    <- collapse_duplicate_genes(expr, symbols)

  # Transformación log2 condicional y universal
  if (needs_log2(expr)) {
    message("  -> Applying log2(x + 0.01) to ", accession)
    expr <- log2(expr + 0.01)
  } else {
    message("  -> ", accession, " already log-scaled; skipping transformation")
  }

  # Ad hoc condition factorization
  level_map <- list(
    GSE76895 = c("ND", "IGT", "T2D"),
    GSE18732 = c("ND", "IGT", "T2D"),
    GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
    GSE27951 = c("NGT", "IGT", "T2D")
  )
 
  exclude_t3cd <- accession == "GSE76895"
  keep <- !is.na(pheno$condition)
  if (exclude_t3cd) keep <- keep & pheno$condition != "T3cD"
  pheno <- pheno[keep, , drop = FALSE]
  expr  <- expr[, pheno$.sample_id, drop = FALSE]
  pheno$condition <- factor(pheno$condition, levels = level_map[[accession]])

  # Per-gene variance (computed before filtering for the global pass)
  gene_var <- apply(expr, 1, var, na.rm = TRUE)

  saveRDS(list(expr = expr, pheno = pheno),
          file.path("data/raw", paste0(accession, "_raw.rds")))
  saveRDS(list(expr = expr, pheno = pheno),
          file.path("data/processed", paste0(accession, "_processed_full.rds")))
  fwrite(pheno,
         file.path("results/qc", paste0(accession, "_pheno.tsv")), sep = "\t")

  summary_df <- pheno |>
    dplyr::count(condition, name = "n") |>
    dplyr::mutate(accession = accession, total_genes = nrow(expr))
  fwrite(summary_df,
         file.path("results/qc", paste0(accession, "_summary.tsv")), sep = "\t")

  invisible(list(expr = expr, pheno = pheno, gene_var = gene_var))
}

raw_list        <- lapply(targets, preprocess_dataset)
names(raw_list) <- targets

# -----------------------------------------------------------------------------
# Global variance filter on shared genes
# -----------------------------------------------------------------------------

# Shared universe
all_gene_sets <- lapply(raw_list, function(x) rownames(x$expr))
common_genes  <- Reduce(intersect, all_gene_sets)
message("Common genes across all datasets: ", length(common_genes))
saveRDS(common_genes, "data/processed/common_genes.rds")
fwrite(data.frame(gene = common_genes), "results/qc/common_genes.tsv", sep = "\t")

# Mean variance per gene over the shared universe
var_matrix    <- do.call(cbind, lapply(raw_list, function(x) x$gene_var[common_genes]))
pooled_var    <- rowMeans(var_matrix, na.rm = TRUE)
names(pooled_var) <- common_genes

global_cutoff <- quantile(pooled_var, 0.5, na.rm = TRUE)
message("Global variance cutoff (Q50): ", round(global_cutoff, 4))

hv_genes <- common_genes[pooled_var >= global_cutoff]
message("High-variance common genes: ", length(hv_genes))
saveRDS(hv_genes, "data/processed/high_variance_genes.rds")
fwrite(data.frame(gene = hv_genes, mean_var = pooled_var[hv_genes]),
       "results/qc/high_variance_genes_global.tsv", sep = "\t")

# Second pass: save processed objects with the global filter
summaries <- lapply(targets, function(acc) {
  obj     <- raw_list[[acc]]
  expr_hv <- obj$expr[hv_genes, , drop = FALSE]
  saveRDS(list(expr = expr_hv, pheno = obj$pheno),
          file.path("data/processed", paste0(acc, "_processed.rds")))
  fwrite(data.frame(gene = rownames(expr_hv)),
         file.path("results/qc", paste0(acc, "_high_variance_genes.tsv")), sep = "\t")
  obj$pheno |>
    dplyr::count(condition, name = "n") |>
    dplyr::mutate(accession = acc, p_genes = nrow(expr_hv))
})

# 1. Unimos todos los resultados en una sola tabla
final_summary <- dplyr::bind_rows(summaries)

# 2. Guardamos el archivo
fwrite(final_summary, "results/qc/dataset_summary_all.tsv", sep = "\t")

# 3. Imprimimos el resultado en la consola
message("\n=================================================")
message("  RESUMEN FINAL: MUESTRAS POR CONDICIÓN Y DATASET")
message("=================================================")
print(as.data.frame(final_summary), row.names = FALSE)
message("=================================================\n")

message("Done.")