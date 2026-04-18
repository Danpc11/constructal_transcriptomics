# Constructal Architecture of the Transcriptome in Type 2 Diabetes

## Theoretical Framework

### Central Hypothesis

Type 2 diabetes (T2D) progression is not a linear degradation of pancreatic islet function but a **topological phase transition**: the co-expression network undergoes a critical reorganization in which the impaired glucose tolerance (IGT) state represents maximum thermodynamic disorder (critical slowing down), and T2D represents a new low-modularity attractor. This hypothesis is tested simultaneously across four metabolically interconnected tissues.

### Constructal Theory Applied to Transcriptomics

Adrian Bejan's constructal law states that any finite-size flow system evolves its architecture over time to facilitate the access of its currents. Applied to gene co-expression networks, this means that healthy transcriptional programs should exhibit configurations that maximize information flow efficiency subject to a metabolic cost constraint. The pipeline operationalizes this through four complementary metrics:

| Metric               | Symbol | Definition                                  | Constructal interpretation               |
|----------------------|--------|---------------------------------------------|------------------------------------------|
| Global efficiency    |  EG    | Mean inverse shortest path length           | Flow accessibility across the network    |
| Communicability      |  Ḡ     | Mean off-diagonal entry of exp(D⁻¹/²WD⁻¹/²) | Propagation efficiency via all paths     |
| Effective resistance |  R̄     | 2·tr(L⁺)/p                                  | Bottleneck cost of information transport |
| Strength entropy     |  Hb    | −Σ pᵢ log pᵢ, pᵢ = kᵢ/Σkⱼ                   | Distributional evenness of connectivity  |

The **Constructal Efficiency Index** (CEI) integrates all four into a single cross-dataset comparable score:

```
CEI = z(EG) + z(Ḡ) − z(R̄) + z(Hb)
```

where z(·) denotes global z-scoring across all datasets and states simultaneously, ensuring cross-tissue comparability.

### Two-Branch Architecture

The pipeline runs two parallel branches that address complementary scientific questions:

**Shared branch (scripts 02–05):** Operates on a fixed subnetwork of up to 800 high-variance genes common to all four datasets. Enables direct cross-tissue comparison of constructal metrics and gene driver scores on identical nodes.

**Full branch (scripts 02b–04b):** Operates on each dataset's complete transcriptome (p = 1,851–21,755 genes). Captures tissue-specific network architecture at full resolution, at the cost of losing cross-dataset node-level comparability.

---

## Datasets: The Metabolic Quartet

| GEO accession | Tissue            | Technology               | n   | States                              |
|---------------|-------------------|--------------------------|-----|-------------------------------------|
| GSE76895      | Pancreatic islets | Affymetrix GPL570        | 83  |(32 ND, 15 IGT, 36 T2D)              |
| GSE18732      | Skeletal muscle   | Affymetrix (Ensembl CDF) | 118 |(47 ND, 26 IGT, 45 T2D)              |
| GSE15653.     | Liver             | Affymetrix GPL570        | 18  |(5 Lean, 4 Obese-noT2D, 9 Obese-T2D) |
| GSE27951      | Adipose tissue    | Affymetrix GPL570        | 33  |(NGT, IGT, T2D)                      |

T3cD samples in GSE76895 are excluded. NGT in GSE27951 is treated as the healthy reference equivalent to ND.

---

## Pipeline Structure

```
01_download_qc_preprocess.R
        │
        ├─ data/processed/{acc}_processed.rds       (shared branch input)
        ├─ data/processed/{acc}_processed_full.rds  (full branch input)
        └─ data/processed/high_variance_genes.rds   (shared subnetwork nodes)
        │
        ├── SHARED BRANCH ──────────────────────────────────────────────────────┐
        │                                                                       │
        ▼                                                                       │
02_networks_modularity_constructal_metrics.R                                    │
        │  Fixed subnetwork (≤800 genes)                                        │
        │  bicor adjacency + WGCNA soft threshold                               │
        │  EG, Ḡ, R̄, Hb, CEI (global z-score)                                   │
        │  Modules: hclust + cutree, maximize n_mod × Q                         │
        └─ results/networks/, results/metrics/                                  │
                │                                                               │
                ▼                                                               │
03_bootstrap_nulls_optimum.R                                                    │
        │  Bootstrap n=100: IC for all 4 metrics                                │
        │  Permutation tests n=1000 (EG+Hb): p_min=0.001                        │
        │  Constructal optimum: W* ∝ (kᵢkⱼ)^(1/α)                               │
        └─ results/bootstrap/, results/nulls/, results/optimum/                 │
                │                                                               │
                ▼                                                               │
04_gene_drivers_and_enrichment.R                                                │
        │  Node metrics: strength, eigencentrality, participation               │
        │  PTI, IRI, TRI (Sherman-Morrison rank-1 update)                       │
        │  KO-support; GO:BP + KEGG enrichment (top-50)                         │
        └─ results/gene_drivers/, results/enrichment/                           │
                │                                                               │
                ├── FULL BRANCH ──────────────────────────────────────────────┐ │
                │                                                             │ │
                ▼                                                             │ │
        02b_full_networks_biology_optimized.R                                 │ │
                │  Complete transcriptome (p = 1k–22k genes)                  │ │
                │  bicor ONE pass + scale-free β selection                    │ │
                │  EG: top-0.1% edges (proportional sparsification)           │ │
                │  Ḡ, R̄: top-3000 hub genes subnetwork                        │ │
                │  Modules: blockwiseModules (WGCNA standard)                 │ │
                └─ results/full_networks/, results/full_metrics/              │ │
                        │                                                     │ │
                        ▼                                                     │ │
                03b_full_bootstrap_and_nulls_optimized.R                      │ │
                        │  Bootstrap n=100 (top-2000 genes): 4 metrics        │ │
                        │  Permutations n=1000 (EG+Hb): p_min=0.001           │ │
                        │  Optimum: Ḡ, R̄ on pre-aligned top-3000 subnet       │ │
                        └─ results/full_bootstrap/, results/full_nulls/,      │ │
                           results/full_optimum/                              │ │
                                │                                             │ │
                                ▼                                             │ │
                        04b_full_gene_drivers_and_enrichment.R                │ │
                                │  PTI/IRI/TRI/KO on top-3000 genes           │ │
                                │  mclapply COW for large W matrices          │ │
                                │  Top-1% edge sparsification for metrics     │ │
                                └─ results/full_gene_drivers/,                │ │
                                   results/full_enrichment/                   │ │
                                                                              │ │
        ◄─────────────────────────────────────────────────────────────────────┘ │
        │                                                                       │
        ◄───────────────────────────────────────────────────────────────────────┘
        │
        ▼
04c_expression_differential.R
        │  limma (correct for RMA-normalized microarrays)
        │  Contrasts: all states vs healthy + IGT→T2D transition
        │  combined_score = rank_pct(PTI+IRI+TRI) × |logFC|
        └─ results/differential_expression/
                │
                ▼
05_cross_dataset_summary.R
        │  Consolidates both branches
        │  Trend tables: EG_change, R̄_change, Ḡ_change, Hb_change
        │  Cross-tissue recurrent driver genes
        └─ results/summary/
```

---

## Mathematical Invariants

These invariants are enforced identically across all scripts and both branches:

**Effective resistance (corrected):**
```
R̄ = 2·tr(L⁺)/p
```
Previous versions used `/(p−1)`, overestimating by a factor of p. Derivation: `R̄ = (2/p(p−1)) · (p−1) · tr(L⁺) = 2·tr(L⁺)/p`.

**Soft-thresholding power (β):** Selected once per state from the bicor matrix by maximizing scale-free fit R² ≥ 0.80 over β ∈ {1,…,20}. Fallback to β = 6 if no β achieves R² ≥ 0.80. Fixed across all bootstrap and permutation iterations — never re-optimized per resample.

**Constructal optimum (first-order approximation):**
```
W*ᵢⱼ ∝ (kᵢ · kⱼ)^(1/α),  subject to Σᵢ<ⱼ (Wᵢⱼ)^α = C
```
where kᵢ are nodal strengths of the healthy reference network and C is its wiring budget. Analytically derived first-order optimum of the Lagrangian `max EG s.t. cost = C`.

**CEI (global z-score):** Always computed across all datasets and states jointly, never within a single dataset. This ensures cross-dataset comparability of sign and magnitude.

---

## Key Methodological Decisions

**bicor vs. Pearson:** Biweight midcorrelation (bicor) is robust to outlier samples common in human tissue microarray data. `maxPOutliers = 0.1` is set consistently across all scripts.

**blockwiseModules in 02b, hclust in 02:** For p ≤ 800 (shared branch), hclust on `1−W` is exact and fast. For p > 3,000 (full branch), TOM-based blockwiseModules is the published WGCNA standard (Langfelder & Horvath 2008). Edge-sparsified fast_greedy was rejected because sparsification can fragment genuine modules.

**EG + Hb for permutation tests, not all 4 metrics:** Gbar (expm, O(p³)) and Rbar (eigen, O(p³)) called 1,000 × 2 × 2 × 4 = 16,000 times would require hours. EG and Hb together detect the same group-level differences under label permutation H0. Gbar and Rbar are reported from observed networks (02b), not the null distribution.

**Rank-percentile for combined_score:** Min-max scaling introduces an edge artifact where the gene with the lowest topological score always gets `combined_score = 0` regardless of |logFC|. Rank normalization to [0,1] is monotonic, has no edge artifacts, and maps naturally to "topological percentile × DE magnitude."

**TRI/KO complexity:** Rank-1 update of the Laplacian pseudoinverse (Sherman-Morrison) reduces complexity from O(p⁴) to O(p²) per gene. Equivalent to full pseudoinverse recomputation to first order in Δw.

**Gbar and Rbar on top-3000 hub genes (02b):** For p = 14,676, `eigen(L, only.values=TRUE)` requires ~1.1T floating-point operations (~10 min). Restricting to the top-3,000 genes by weighted degree (strength) reduces this to ~27B ops (~3 seconds) with <5% empirical error, as spectral properties are dominated by hub connectivity.

---

## Computational Requirements

| Script |   RAM  | CPU (40 workers) | Bottleneck                      |
|--------|--------|------------------|---------------------------------|
| 01     | ~4 GB  | ~20 min          | GEO download                    |
| 02     | ~2 GB  | ~5 min           | bicor (p≤800)                   |
| 02b    | ~32 GB | ~30–60 min       | bicor (p=14k), blockwiseModules |
| 03     | ~4 GB  | ~10 min          | 100 bootstrap × 4 metrics       |
| 03b    | ~16 GB | ~20 min          | 1000 permutations (EG+Hb)       |
| 04     | ~8 GB  | ~5 min           | TRI rank-1 (p≤800)              |
| 04b    | ~32 GB | ~15 min          | TRI rank-1 (p≤1500)             |
| 04c    | ~4 GB  | ~5 min           | limma + enrichGO                |
| 05     | ~2 GB  | <1 min           | File consolidation              |

### Required R packages

```r
# CRAN
install.packages(c(
  "WGCNA", "data.table", "dplyr", "tidyr", "igraph",
  "expm", "Matrix", "parallel", "foreach", "doParallel",
  "stringr", "tibble"
))

# Recommended (graceful fallback if absent)
install.packages(c("RSpectra", "doRNG"))

# Bioconductor
BiocManager::install(c(
  "limma", "GEOquery", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "enrichplot",
  "hthgu133a.db"
))
```

---

## Execution Order

```bash
# 1. Download, QC, normalize (required first)
Rscript 01_download_qc_preprocess.R

# 2a. Shared branch (sequential)
Rscript 02_networks_modularity_constructal_metrics.R
Rscript 03_bootstrap_nulls_optimum.R
Rscript 04_gene_drivers_and_enrichment.R

# 2b. Full branch (HPC; can run in parallel with 2a)
Rscript 02b_full_networks_biology_optimized.R
Rscript 03b_full_bootstrap_and_nulls_optimized.R
Rscript 04b_full_gene_drivers_and_enrichment.R

# 3. Differential expression + integration (requires both branches)
Rscript 04c_expression_differential.R

# 4. Cross-dataset summary (requires all previous steps)
Rscript 05_cross_dataset_summary.R
```

Individual datasets can be targeted: `Rscript 02_networks.R GSE76895 GSE18732`

---

## Output Files

```
results/
├── qc/                         # Sample counts, phenotype tables, gene lists
├── metrics/                    # Shared: EG, Ḡ, R̄, Hb, CEI per state
├── full_metrics/               # Full: same metrics at full resolution
├── networks/                   # Shared: adjacency matrices W, β
├── full_networks/              # Full: adjacency matrices W, β
├── modules/                    # Shared: module membership per gene
├── bootstrap/                  # Shared: bootstrap metric distributions
├── nulls/                      # Shared: permutation p-values (EG, Hb, CEI)
├── optimum/                    # Shared: W*, ΔE, ΔR, ΔG, ΔH
├── full_bootstrap/             # Full branch equivalents
├── full_nulls/
├── full_optimum/
├── gene_drivers/               # PTI, IRI, TRI, KO, class, GO/KEGG
├── full_gene_drivers/
├── enrichment/
├── full_enrichment/
├── differential_expression/    # limma results + combined_score integration
└── summary/                    # Cross-dataset consolidated tables
```

### Key output files

| File | Content |
|------|---------|
| `results/metrics/all_constructal_metrics.tsv`           | Main constructal metrics (shared branch) |
| `results/full_metrics/all_full_constructal_metrics.tsv` | Full-resolution metrics                  |
| `results/summary/cross_dataset_trends.tsv`              | EG_change, R̄_change per tissue/branch    |
| `results/summary/all_deviation_from_optimum.tsv`        | ΔE, ΔR, ΔG, ΔH relative to W*            |
| `results/gene_drivers/{acc}_gene_driver_scores.tsv`     | Per-gene PTI/IRI/TRI/class               |
| `results/differential_expression/{acc}_full_integrated_drivers_DE.tsv` | combined_score            |
| `results/summary/all_integrated_drivers_DE.tsv`         | Cross-tissue integrated ranking          |

---

## Gene Driver Score Definitions

| Score | Full name              | Definition                  | High value means |
|-------|------------------------|-----------------------------|------------------|
| PTI | Phase Transition Index   | Node topology change ND→IGT | Gene rewires most at the critical transition |
| IRI | Irreversibility Index    | Node topology change ND→T2D | Gene is most altered in established disease  |
| TRI | Topological Rescue Index | ΔEG + ΔHb if gene's edges restored to healthy | Restoring this gene improves network most.  |
| KO  | KO-support               | ΔΔEG if gene is knocked out | Gene sustains the diseased network topology  |
| combined_score | Integrated driver | rank_pct(PTI+IRI+TRI) × \|logFC\| | Topologically important AND differentially expressed |

Gene classes (mutually exclusive, evaluated in priority order):

| Class                    | Criteria            | Biological interpretation                   |
|--------------------------|---------------------|---------------------------------------------|
| `transition_rescue`      | high PTI + high TRI | Early driver with rescue potential          |
| `irreversibility_rescue` | high IRI + high TRI | Late driver with rescue potential           |
| `irreversibility_lock`   | high IRI + high KO  | Sustains the diseased attractor             |
| `transition_driver`      | high PTI            | Reorganization driver without rescue effect |
| `other`                  | none of the above   | Background topology                         |

---

## Reproducibility

- All parallel scripts: `set.seed(1234)` with `RNGkind("L'Ecuyer-CMRG")`
- doRNG (if installed) provides per-iteration reproducibility in foreach loops
- mclapply uses L'Ecuyer-CMRG streams automatically via `set.seed()` before forking
- β values are saved in network `.rds` files and reused identically in all downstream scripts
- All scripts accept dataset targets as command-line arguments for partial re-runs

---

## Author 

- Daniel Pérez Calixto -

Instituto Nacional de Médicina Genómica

Contact info: dperez@inmegen.gob.mx
