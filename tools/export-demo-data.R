#!/usr/bin/env Rscript
# Export precomputed demo data for the interactive UI (NCI-60 + Golub).
# Run via Docker: docker-compose run --rm evaluomeR bash tools/docker-export.sh
# Expected runtime: ~30-50 min (both datasets, bs=100)

suppressPackageStartupMessages({
  library(evaluomeR)
  library(SummarizedExperiment)
  library(jsonlite)
  library(mclust)
  library(cluster)
})

BS    <- 100L
SEED  <- 13606L
TOP_N_CHART_METRICS <- 8L

as_metric_df <- function(data_obj) {
  if (inherits(data_obj, "SummarizedExperiment")) {
    as.data.frame(assay(data_obj))
  } else {
    as.data.frame(data_obj)
  }
}

run_export <- function(dataset_name,
                       data_obj,
                       k_min,
                       k_max,
                       bs,
                       seed,
                       cancer_type_fn,
                       out_file,
                       sample_id_fn = NULL) {
  raw_df <- as_metric_df(data_obj)
  message("\n=== Exporting ", dataset_name, " ===")
  message("Samples: ", nrow(raw_df), " · metrics: ", ncol(raw_df) - 1)

  raw_descriptions <- as.character(raw_df$Description)
  cancer_types <- vapply(raw_descriptions, cancer_type_fn, character(1))
  unique_types <- sort(unique(cancer_types))
  message("Labels: ", paste(unique_types, collapse = ", "))

  if (is.null(sample_id_fn)) {
    sample_ids <- raw_descriptions
  } else {
    sample_ids <- sample_id_fn(cancer_types, raw_descriptions)
  }

  metric_cols <- setdiff(colnames(raw_df), c("Description", "Class"))
  metric_df   <- raw_df[, metric_cols, drop = FALSE]

  # Coerce to numeric (golub may ship factor columns)
  metric_df <- as.data.frame(lapply(metric_df, function(x) suppressWarnings(as.numeric(as.character(x)))))
  metric_cols <- names(metric_df)[vapply(metric_df, is.numeric, logical(1))]
  metric_df   <- metric_df[, metric_cols, drop = FALSE]

  # Data frame passed to evaluomeR (metrics only, no Class)
  analysis_df <- cbind(Description = raw_df$Description, metric_df)

  col_vars <- vapply(metric_df, var, numeric(1), na.rm = TRUE)
  chart_metrics <- names(sort(col_vars, decreasing = TRUE))[seq_len(
    min(TOP_N_CHART_METRICS, length(col_vars))
  )]

  samples <- lapply(seq_len(nrow(raw_df)), function(i) {
    m_vals <- as.list(unname(as.numeric(metric_df[i, ])))
    names(m_vals) <- metric_cols
    list(
      id         = sample_ids[i],
      cancerType = unname(cancer_types[i]),
      metrics    = m_vals
    )
  })

  message("stabilityRange k=", k_min, ":", k_max, " bs=", bs, " all_metrics=TRUE ...")
  stab_range <- stabilityRange(
    analysis_df,
    k.range     = c(k_min, k_max),
    bs          = bs,
    getImages   = FALSE,
    all_metrics = TRUE,
    seed        = seed
  )
  stab_std <- standardizeStabilityData(stab_range)

  stab_list <- lapply(k_min:k_max, function(k) {
    col <- paste0("k_", k)
    setNames(as.list(stab_std[, col]), rownames(stab_std))
  })
  names(stab_list) <- as.character(k_min:k_max)

  message("qualityRange (silhouette) all_metrics=TRUE ...")
  qual_range <- qualityRange(
    analysis_df,
    k.range     = c(k_min, k_max),
    getImages   = FALSE,
    all_metrics = TRUE,
    seed        = seed
  )
  qual_std <- standardizeQualityData(qual_range)

  sil_list <- lapply(k_min:k_max, function(k) {
    col <- paste0("k_", k)
    setNames(as.list(qual_std[, col]), rownames(qual_std))
  })
  names(sil_list) <- as.character(k_min:k_max)

  message("qualityRange (CH) all_metrics=TRUE ...")
  ch_data <- qualityRange(
    analysis_df,
    k.range       = c(k_min, k_max),
    getImages     = FALSE,
    all_metrics   = TRUE,
    seed          = seed,
    quality_index = "ch"
  )
  ch_list <- lapply(k_min:k_max, function(k) {
    df <- ch_data[[paste0("k=", k)]]
    setNames(as.list(df$CH), df$Metric)
  })
  names(ch_list) <- as.character(k_min:k_max)

  message("getOptimalKValue ...")
  opt_df <- getOptimalKValue(stab_range, qual_range)
  opt_per_metric <- setNames(
    as.list(as.integer(opt_df$Global_optimal_k)),
    opt_df$Metric
  )
  global_opt_k <- as.integer(opt_df$Global_optimal_k[1])
  message("Global optimal k = ", global_opt_k)

  opt_detail <- lapply(seq_len(nrow(opt_df)), function(i) {
    list(
      metric         = as.character(opt_df$Metric[i]),
      stabilityMaxK  = as.integer(opt_df$Stability_max_k[i]),
      qualityMaxK    = as.integer(opt_df$Quality_max_k[i]),
      globalOptimalK = as.integer(opt_df$Global_optimal_k[i])
    )
  })

  # Aggregated all_metrics index (single row from evaluomeR)
  index_metrics <- intersect("all_metrics", rownames(stab_std))
  if (length(index_metrics) == 0) index_metrics <- rownames(stab_std)

  clusters_list <- lapply(k_min:k_max, function(k) {
    set.seed(seed)
    km <- kmeans(metric_df, centers = k, nstart = 25)
    as.integer(km$cluster)
  })
  names(clusters_list) <- as.character(k_min:k_max)

  message("ARI vs ground truth per k ...")
  true_int <- as.integer(as.factor(cancer_types))
  ari_by_k <- setNames(
    lapply(k_min:k_max, function(k) {
      round(mclust::adjustedRandIndex(true_int, clusters_list[[as.character(k)]]), 6)
    }),
    as.character(k_min:k_max)
  )

  message("Stability SD per k (inter-cluster Jaccard SD) ...")
  stab_sd_by_k <- setNames(
    lapply(k_min:k_max, function(k) {
      raw <- clusterbootWrapper(
        data = metric_df,
        B = bs,
        bootmethod = "boot",
        cbi = "kmeans",
        krange = k,
        gold_standard = NULL,
        seed = seed
      )
      round(sd(raw$bootmean, na.rm = TRUE), 6)
    }),
    as.character(k_min:k_max)
  )

  message("Per-sample silhouette per k ...")
  scaled_df <- scale(metric_df)
  dist_mat <- dist(scaled_df)
  sil_samples_by_k <- setNames(
    lapply(k_min:k_max, function(k) {
      sil <- silhouette(clusters_list[[as.character(k)]], dist_mat)
      round(as.numeric(sil[, "sil_width"]), 6)
    }),
    as.character(k_min:k_max)
  )

  k_summary <- lapply(k_min:k_max, function(k) {
    k_col <- paste0("k_", k)
    avg_stab <- as.numeric(unlist(stab_std[index_metrics, k_col]))
    avg_sil  <- as.numeric(unlist(qual_std[index_metrics, k_col]))
    composite <- avg_stab * avg_sil
    list(
      avgStability  = round(avg_stab, 6),
      avgSilhouette = round(avg_sil, 6),
      composite     = round(composite, 6),
      ari           = ari_by_k[[as.character(k)]]
    )
  })
  names(k_summary) <- as.character(k_min:k_max)

  message("PCA (3 components) ...")
  pca_result    <- prcomp(metric_df, center = TRUE, scale. = TRUE)
  n_pcs         <- min(3L, ncol(pca_result$x))
  pc_coords     <- as.data.frame(pca_result$x[, seq_len(n_pcs), drop = FALSE])
  var_explained <- summary(pca_result)$importance[2, seq_len(n_pcs)]
  var_msg <- paste(
    vapply(seq_len(n_pcs), function(i) {
      paste0("PC", i, "=", round(var_explained[i] * 100, 1), "%")
    }, character(1)),
    collapse = ", "
  )
  message("PCA variance: ", var_msg)

  pca_export <- list(
    coords = lapply(seq_len(nrow(pc_coords)), function(i) {
      coord <- list(
        pc1 = round(pc_coords$PC1[i], 6),
        pc2 = round(pc_coords$PC2[i], 6)
      )
      if (n_pcs >= 3L) {
        coord$pc3 <- round(pc_coords$PC3[i], 6)
      }
      coord
    }),
    varExplained = as.numeric(var_explained)
  )

  demo <- list(
    meta = list(
      dataset      = dataset_name,
      nSamples     = nrow(raw_df),
      nCancerTypes = length(unique_types),
      cancerTypes  = as.list(unique_types),
      metrics        = metric_cols,
      chartMetrics   = as.list(index_metrics),
      previewMetrics = as.list(chart_metrics[seq_len(min(3L, length(chart_metrics)))]),
      allMetrics     = TRUE,
      kRange       = c(k_min, k_max),
      bs           = bs
    ),
    samples           = samples,
    stability         = stab_list,
    quality = list(
      silhouette = sil_list,
      ch         = ch_list
    ),
    optimalK          = global_opt_k,
    optimalKPerMetric = opt_per_metric,
    optimalKDetail    = opt_detail,
    kSummary          = k_summary,
    stabilitySD       = stab_sd_by_k,
    sampleSilhouette  = sil_samples_by_k,
    clusters          = clusters_list,
    pca               = pca_export
  )

  out_dir <- file.path("ui", "src", "data")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(demo, out_file, auto_unbox = TRUE, digits = 6, pretty = FALSE)
  message("Written ", out_file)
}

# ── NCI-60 ───────────────────────────────────────────────────────────────────
data("nci60_k8")

nci_cancer_type <- function(d) toupper(trimws(d))
nci_sample_ids <- function(cancer_types, descriptions) {
  unique_types <- sort(unique(cancer_types))
  type_counter <- setNames(integer(length(unique_types)), unique_types)
  ids <- character(length(cancer_types))
  for (i in seq_along(cancer_types)) {
    ct <- cancer_types[i]
    type_counter[ct] <- type_counter[ct] + 1L
    ids[i] <- paste0(ct, "_", type_counter[ct])
  }
  ids
}

run_export(
  dataset_name   = "nci60_k8",
  data_obj       = nci60_k8,
  k_min          = 3L,
  k_max          = 8L,
  bs             = BS,
  seed           = SEED,
  cancer_type_fn = nci_cancer_type,
  out_file       = file.path("ui", "src", "data", "nci60.json"),
  sample_id_fn   = nci_sample_ids
)

# ── Golub leukemia ───────────────────────────────────────────────────────────
data("golub")

golub_cancer_type <- function(d) {
  s <- sub("[0-9]+$", "", toupper(trimws(d)))
  if (s == "B") return("B-ALL")
  if (s == "T") return("T-ALL")
  s
}

run_export(
  dataset_name   = "golub",
  data_obj       = golub,
  k_min          = 2L,
  k_max          = 6L,
  bs             = BS,
  seed           = SEED,
  cancer_type_fn = golub_cancer_type,
  out_file       = file.path("ui", "src", "data", "golub.json"),
  sample_id_fn   = NULL
)
