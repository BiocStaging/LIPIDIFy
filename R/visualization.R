#' Visualize Raw Data
#'
#' Produces a simple per-sample boxplot, density, or histogram of raw
#' lipidomics intensities.
#'
#' @param data_list List as returned by \code{load_lipidomics_data} or
#'   \code{load_lipidomics_data_from_df}, must contain \code{$numeric_data}.
#' @param plot_type One of \code{"boxplot"}, \code{"density"}, or
#'   \code{"histogram"}.
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' dl <- load_lipidomics_data_from_df(generate_example_data())
#' p <- visualize_raw_data(dl, "boxplot")
#' print(p)
visualize_raw_data <- function(data_list, plot_type = "boxplot") {
  if (!is.list(data_list) || is.null(data_list$numeric_data)) {
    stop(
      "`data_list` must be a list containing a `numeric_data` element ",
      "(e.g. the output of load_lipidomics_data() or load_lipidomics_data_from_df())."
    )
  }
  numeric_data <- data_list$numeric_data
  if (is.matrix(numeric_data)) {
    numeric_data <- as.data.frame(numeric_data, check.names = FALSE)
  }

  to_long <- function(d) {
    long <- tidyr::pivot_longer(d, dplyr::everything(),
      names_to = "Sample", values_to = "Intensity"
    )
    long$Intensity <- as.numeric(long$Intensity)
    long[!is.na(long$Intensity), ]
  }

  long <- to_long(numeric_data)

  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(long, ggplot2::aes(x = Sample, y = Intensity)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Raw Data Distribution", y = "Intensity")
  } else if (plot_type == "density") {
    p <- ggplot2::ggplot(long, ggplot2::aes(x = Intensity, color = Sample)) +
      ggplot2::geom_density(alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "Raw Data Density Distribution", x = "Intensity")
  } else if (plot_type == "histogram") {
    p <- ggplot2::ggplot(long, ggplot2::aes(x = Intensity)) +
      ggplot2::geom_histogram(bins = 50, alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Raw Data Histogram", x = "Intensity", y = "Count")
  } else {
    stop("plot_type must be one of 'boxplot', 'density', or 'histogram'.")
  }

  p
}


#' Visualize Raw or Normalized Data with Sample/Lipid Toggle
#'
#' Extended visualization that can show data from the sample perspective
#' (one boxplot/violin/density per sample) or the lipid perspective (top
#' variable lipids).
#'
#' @param data_list List with \code{$numeric_data} (samples \eqn{\times}
#'   lipids matrix).
#' @param plot_type One of \code{"boxplot"}, \code{"violin"}, \code{"density"},
#'   or \code{"histogram"}.
#' @param view_mode Either \code{"sample"} or \code{"lipid"}.
#' @param top_n Integer or \code{NULL}. In lipid mode, the number of most
#'   variable lipids to display (defaults to 30 when \code{NULL}). In sample
#'   mode, the number of most variable samples to display; \code{NULL} (the
#'   default) shows all samples.
#' @param metadata Optional metadata data frame (same row order as
#'   \code{numeric_data}) used to colour samples by group when
#'   \code{group_column} is provided.
#' @param group_column Name of the group column in \code{metadata}.
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' dl <- load_lipidomics_data_from_df(generate_example_data())
#' p <- visualize_raw_data_improved(dl, "boxplot", "sample")
#' print(p)
visualize_raw_data_improved <- function(data_list,
                                        plot_type = "boxplot",
                                        view_mode = "sample",
                                        top_n = NULL,
                                        metadata = NULL,
                                        group_column = "Sample Group") {
  numeric_data <- data_list$numeric_data
  if (is.null(numeric_data)) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = "No numeric_data found"
      ) +
      ggplot2::theme_void())
  }

  # Build a data frame with a Sample column from row names
  df <- as.data.frame(numeric_data, check.names = FALSE)
  if (is.null(rownames(df))) rownames(df) <- paste0("S", seq_len(nrow(df)))
  df$Sample <- rownames(df)

  # Attach group information if metadata is available
  group_col_present <- FALSE
  if (!is.null(metadata) && group_column %in% colnames(metadata)) {
    group_col_present <- TRUE
    df$Group <- as.character(metadata[[group_column]])
  }

  # Helper: pivot to long format and coerce Intensity
  pivot_long_safe <- function(d, id_vars) {
    long <- tidyr::pivot_longer(d,
      cols = -tidyr::all_of(id_vars),
      names_to = "Lipid", values_to = "Intensity"
    )
    long$Intensity <- as.numeric(long$Intensity)
    long[!is.na(long$Intensity), ]
  }

  plot_type <- match.arg(plot_type, c("boxplot", "violin", "density", "histogram"))

  # -------- Sample view --------------------------------------------------
  if (identical(view_mode, "sample")) {
    id_vars <- if (group_col_present) c("Sample", "Group") else "Sample"

    # Optionally restrict to the top_n most variable samples (variance taken
    # across each sample's lipid intensities). NULL/NA keeps all samples.
    if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
      feature_cols <- setdiff(colnames(df), c("Sample", "Group"))
      if (length(feature_cols) > 0 && nrow(df) > top_n) {
        sample_vars <- apply(df[, feature_cols, drop = FALSE], 1,
          stats::var,
          na.rm = TRUE
        )
        sample_vars[is.na(sample_vars)] <- 0
        keep <- order(sample_vars, decreasing = TRUE)[seq_len(top_n)]
        df <- df[sort(keep), , drop = FALSE]
      }
    }

    long <- pivot_long_safe(df, id_vars)

    # Sort samples by group if group information is available
    if (group_col_present) {
      sample_order <- unique(long$Sample[order(long$Group, long$Sample)])
      long$Sample <- factor(long$Sample, levels = sample_order)
    }

    fill_aes <- if (group_col_present) {
      ggplot2::aes(
        x = Sample, y = Intensity,
        fill = Group
      )
    } else {
      ggplot2::aes(x = Sample, y = Intensity)
    }

    color_aes <- if (group_col_present) {
      ggplot2::aes(x = Intensity, color = Group)
    } else {
      ggplot2::aes(x = Intensity, color = Sample)
    }

    if (plot_type == "boxplot") {
      p <- ggplot2::ggplot(long, fill_aes) +
        ggplot2::geom_boxplot(outlier.alpha = 0.3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) +
        ggplot2::labs(title = "Sample-wise Distribution", x = "Sample", y = "Intensity")
    } else if (plot_type == "violin") {
      p <- ggplot2::ggplot(long, fill_aes) +
        ggplot2::geom_violin(trim = FALSE) +
        ggplot2::geom_boxplot(
          width = 0.12, outlier.alpha = 0.2,
          fill = "white", color = "grey40"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) +
        ggplot2::labs(
          title = "Sample-wise Distribution (Violin)",
          x = "Sample", y = "Intensity"
        )
    } else if (plot_type == "density") {
      p <- ggplot2::ggplot(long, color_aes) +
        ggplot2::geom_density(alpha = 0.7) +
        ggplot2::theme_minimal() +
        (if (!group_col_present) ggplot2::theme(legend.position = "none") else NULL) +
        ggplot2::labs(title = "Sample-wise Density", x = "Intensity", y = "Density")
    } else { # histogram
      p <- ggplot2::ggplot(long, ggplot2::aes(x = Intensity)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Sample-wise Histogram",
          x = "Intensity", y = "Count"
        )
    }

    # -------- Lipid view ---------------------------------------------------
  } else if (identical(view_mode, "lipid")) {
    if (is.null(top_n)) top_n <- 30

    feature_cols <- setdiff(colnames(df), c("Sample", "Group"))
    vars <- apply(df[, feature_cols, drop = FALSE], 2,
      stats::var,
      na.rm = TRUE
    )
    vars[is.na(vars)] <- 0
    top_lipids <- names(sort(vars, decreasing = TRUE))[
      seq_len(min(top_n, length(vars)))
    ]

    if (length(top_lipids) == 0) {
      return(ggplot2::ggplot() +
        ggplot2::annotate("text",
          x = 0.5, y = 0.5,
          label = "No lipid columns found"
        ) +
        ggplot2::theme_void())
    }

    id_vars <- if (group_col_present) c("Sample", "Group") else "Sample"
    df_sub <- df[, c(id_vars, top_lipids), drop = FALSE]
    long <- pivot_long_safe(df_sub, id_vars)

    fill_aes <- if (group_col_present) {
      ggplot2::aes(
        x = Lipid, y = Intensity,
        fill = Group
      )
    } else {
      ggplot2::aes(x = Lipid, y = Intensity)
    }

    if (plot_type == "boxplot") {
      p <- ggplot2::ggplot(long, fill_aes) +
        ggplot2::geom_boxplot(outlier.alpha = 0.3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)) +
        ggplot2::labs(
          title = paste("Top", length(top_lipids), "Variable Lipids"),
          x = "Lipid", y = "Intensity"
        )
    } else if (plot_type == "violin") {
      p <- ggplot2::ggplot(long, fill_aes) +
        ggplot2::geom_violin(trim = FALSE) +
        ggplot2::geom_boxplot(
          width = 0.12, outlier.alpha = 0.2,
          fill = "white", color = "grey40"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)) +
        ggplot2::labs(
          title = paste("Top", length(top_lipids), "Variable Lipids (Violin)"),
          x = "Lipid", y = "Intensity"
        )
    } else if (plot_type == "density") {
      color_aes2 <- if (group_col_present) {
        ggplot2::aes(x = Intensity, color = Group)
      } else {
        ggplot2::aes(x = Intensity, color = Lipid)
      }
      p <- ggplot2::ggplot(long, color_aes2) +
        ggplot2::geom_density(alpha = 0.7) +
        ggplot2::theme_minimal() +
        (if (!group_col_present) ggplot2::theme(legend.position = "none") else NULL) +
        ggplot2::labs(
          title = paste("Density - Top", length(top_lipids), "Variable Lipids"),
          x = "Intensity", y = "Density"
        )
    } else { # histogram
      p <- ggplot2::ggplot(long, ggplot2::aes(x = Intensity)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.7, fill = "darkgreen") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Histogram - Top", length(top_lipids), "Variable Lipids"),
          x = "Intensity", y = "Count"
        )
    }
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = paste("Unknown view mode:", view_mode)
      ) +
      ggplot2::theme_void()
  }

  p
}


#' Quick QC Plot for a Normalized Data Matrix
#'
#' @param data_matrix Numeric matrix or data frame (samples as rows,
#'   lipids as columns).
#' @param title Plot title string.
#' @param metadata Optional metadata data frame (same row order) for group
#'   colouring.
#' @param group_column Name of the group column in \code{metadata}.
#' @param plot_type One of \code{"boxplot"}, \code{"violin"}, or
#'   \code{"density"}. Defaults to \code{"boxplot"}.
#' @param view_mode Either \code{"sample"} (one box/violin per sample; the
#'   default, preserving previous behaviour) or \code{"lipid"} (one box/violin
#'   per lipid). Ignored for \code{"density"}, which always shows the intensity
#'   distribution.
#' @param top_n Integer or \code{NULL}. In \code{"sample"} view, restrict the
#'   plot to the \code{top_n} most variable samples (\code{NULL}, the default,
#'   shows all). In \code{"lipid"} view, restrict to the \code{top_n} most
#'   variable lipids (\code{NULL} defaults to 30).
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' m <- matrix(rlnorm(60, 8, 1), nrow = 6, ncol = 10)
#' p <- create_pipeline_plot(m, title = "Test Pipeline")
#' print(p)
create_pipeline_plot <- function(data_matrix, title = "Normalization Pipeline",
                                 metadata = NULL, group_column = "Sample Group",
                                 plot_type = "boxplot",
                                 view_mode = "sample", top_n = NULL) {
  plot_type <- match.arg(plot_type, c("boxplot", "violin", "density"))
  view_mode <- match.arg(view_mode, c("sample", "lipid"))
  .validate_numeric_matrix(data_matrix, "data_matrix")
  # group_column is optional here by design: if absent from metadata, the
  # plot is simply drawn without group colouring (group_col_present below),
  # so only the data.frame-ness of metadata itself is validated.
  if (!is.null(metadata)) .validate_metadata(metadata)

  df <- as.data.frame(data_matrix, check.names = FALSE)

  if (is.null(rownames(df)) ||
    all(rownames(df) == as.character(seq_len(nrow(df))))) {
    rownames(df) <- paste0("Sample_", seq_len(nrow(df)))
  }
  df$Sample <- rownames(df)

  group_col_present <- FALSE
  if (!is.null(metadata) && group_column %in% colnames(metadata)) {
    group_col_present <- TRUE
    df$Group <- as.character(metadata[[group_column]])
  }

  id_vars <- if (group_col_present) c("Sample", "Group") else "Sample"

  # ---- Optional top-N selection (most variable samples or lipids) ----------
  if (identical(view_mode, "lipid")) {
    tn <- if (is.null(top_n) || is.na(top_n)) 30L else as.integer(top_n)
    feature_cols <- setdiff(colnames(df), c("Sample", "Group"))
    if (length(feature_cols) > tn) {
      vars <- apply(df[, feature_cols, drop = FALSE], 2, stats::var, na.rm = TRUE)
      vars[is.na(vars)] <- 0
      keep_lipids <- names(sort(vars, decreasing = TRUE))[seq_len(tn)]
      df <- df[, c(id_vars, keep_lipids), drop = FALSE]
    }
  } else if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
    feature_cols <- setdiff(colnames(df), c("Sample", "Group"))
    if (length(feature_cols) > 0 && nrow(df) > top_n) {
      sample_vars <- apply(df[, feature_cols, drop = FALSE], 1, stats::var, na.rm = TRUE)
      sample_vars[is.na(sample_vars)] <- 0
      keep <- order(sample_vars, decreasing = TRUE)[seq_len(top_n)]
      df <- df[sort(keep), , drop = FALSE]
    }
  }

  long <- tidyr::pivot_longer(df,
    cols = -tidyr::all_of(id_vars),
    names_to = "Lipid", values_to = "Intensity"
  )
  long$Intensity <- as.numeric(long$Intensity)
  long <- long[!is.na(long$Intensity), ]

  if (group_col_present) {
    sample_order <- unique(long$Sample[order(long$Group, long$Sample)])
    long$Sample <- factor(long$Sample, levels = sample_order)
  }

  # x axis is samples or lipids depending on the view (density is unaffected)
  x_col <- if (identical(view_mode, "lipid")) "Lipid" else "Sample"
  x_lab <- x_col

  fill_aes <- if (group_col_present) {
    ggplot2::aes(x = .data[[x_col]], y = Intensity, fill = Group)
  } else {
    ggplot2::aes(x = .data[[x_col]], y = Intensity)
  }
  color_aes <- if (group_col_present) {
    ggplot2::aes(x = Intensity, color = Group)
  } else {
    ggplot2::aes(x = Intensity, color = Sample)
  }

  base_theme <- list(
    ggplot2::theme_minimal(),
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  )

  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(long, fill_aes) +
      ggplot2::geom_boxplot(outlier.alpha = 0.3) +
      base_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)) +
      ggplot2::labs(title = title, x = x_lab, y = "Intensity")
  } else if (plot_type == "violin") {
    p <- ggplot2::ggplot(long, fill_aes) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
      ggplot2::geom_boxplot(
        width = 0.1, outlier.alpha = 0.2,
        fill = "white", color = "grey40"
      ) +
      base_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)) +
      ggplot2::labs(title = title, x = x_lab, y = "Intensity")
  } else {
    p <- ggplot2::ggplot(long, color_aes) +
      ggplot2::geom_density(alpha = 0.5) +
      base_theme +
      (if (!group_col_present) ggplot2::theme(legend.position = "none") else NULL) +
      ggplot2::labs(title = title, x = "Intensity", y = "Density")
  }

  p
}


#' Create a Robust Heatmap of Top Variable Features
#'
#' @param data_matrix Numeric matrix (features as rows, samples as columns).
#' @param metadata Metadata data frame (samples as rows).
#' @param group_column Name of the group column in \code{metadata}.
#' @param top_n Maximum number of features to display.
#' @param classification_data Optional classification data frame for row
#'   annotation (must have a \code{Lipid} column).
#' @param title Heatmap title.
#' @return A \code{pheatmap} object, or a \code{ggplot2} error plot.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' create_heatmap_robust(t(norm), d$metadata, "Sample Group", top_n = 10)
create_heatmap_robust <- function(data_matrix, metadata,
                                  group_column = "Sample Group",
                                  top_n = 50,
                                  classification_data = NULL,
                                  title = "Heatmap") {
  tryCatch(
    {
      if (is.null(data_matrix) || nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
        stop("Data matrix is empty or NULL.")
      }
      .validate_metadata(metadata)

      # Orient with features as rows / samples as columns, resolved by sample
      # identity rather than by comparing dimension sizes.
      data_matrix <- .orient_matrix(data_matrix, metadata, want = "features_rows")

      # Resolve sample names
      sample_names <- if ("Sample Name" %in% colnames(metadata)) {
        metadata$`Sample Name`
      } else if (!is.null(rownames(metadata))) {
        rownames(metadata)
      } else {
        paste0("Sample_", seq_len(nrow(metadata)))
      }

      # Align samples
      common_samples <- intersect(colnames(data_matrix), sample_names)
      if (length(common_samples) == 0) {
        n_s <- min(ncol(data_matrix), nrow(metadata))
        data_matrix <- data_matrix[, seq_len(n_s), drop = FALSE]
        metadata <- metadata[seq_len(n_s), , drop = FALSE]
        colnames(data_matrix) <- sample_names[seq_len(n_s)]
      } else {
        data_matrix <- data_matrix[, common_samples, drop = FALSE]
        metadata <- metadata[sample_names %in% common_samples, , drop = FALSE]
      }

      # Select top variable features
      if (nrow(data_matrix) > top_n) {
        feature_vars <- apply(data_matrix, 1, stats::var, na.rm = TRUE)
        top_features <- names(sort(feature_vars, decreasing = TRUE)[seq_len(top_n)])
        data_matrix <- data_matrix[top_features, , drop = FALSE]
      }

      scaled_data <- t(scale(t(data_matrix)))
      scaled_data[is.na(scaled_data)] <- 0

      annotation_col <- if (group_column %in% colnames(metadata)) {
        data.frame(Group = metadata[[group_column]], row.names = colnames(data_matrix))
      } else {
        NA
      }

      pheatmap::pheatmap(
        scaled_data,
        annotation_col = annotation_col,
        main           = title,
        cluster_rows   = TRUE,
        cluster_cols   = TRUE,
        show_rownames  = nrow(scaled_data) <= 50,
        show_colnames  = TRUE,
        fontsize_row   = 8,
        fontsize_col   = 8,
        silent         = TRUE
      )
    },
    error = function(e) {
      ggplot2::ggplot() +
        ggplot2::annotate("text",
          x = 0.5, y = 0.5,
          label = paste("Heatmap failed:\n", e$message), size = 4
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Heatmap Error")
    }
  )
}


#' Align Samples Between Data Matrix and Metadata
#'
#' @param data_matrix Numeric matrix (features as rows, samples as columns).
#' @param metadata Metadata data frame.
#' @return Named list with aligned \code{data_matrix} and \code{metadata}.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' aligned <- fix_sample_alignment(t(d$numeric_data), d$metadata)
#' names(aligned)
fix_sample_alignment <- function(data_matrix, metadata) {
  .validate_metadata(metadata)
  data_samples <- colnames(data_matrix)

  metadata_samples <- if ("Sample Name" %in% colnames(metadata)) {
    metadata$`Sample Name`
  } else if (!is.null(rownames(metadata))) {
    rownames(metadata)
  } else {
    paste0("Sample_", seq_len(nrow(metadata)))
  }

  common_samples <- intersect(data_samples, metadata_samples)

  if (length(common_samples) > 0) {
    data_matrix <- data_matrix[, common_samples, drop = FALSE]
    if ("Sample Name" %in% colnames(metadata)) {
      metadata <- metadata[metadata$`Sample Name` %in% common_samples, , drop = FALSE]
      rownames(metadata) <- metadata$`Sample Name`
    } else {
      metadata <- metadata[rownames(metadata) %in% common_samples, , drop = FALSE]
    }
    metadata <- metadata[colnames(data_matrix), , drop = FALSE]
  } else {
    n_s <- min(ncol(data_matrix), nrow(metadata))
    data_matrix <- data_matrix[, seq_len(n_s), drop = FALSE]
    metadata <- metadata[seq_len(n_s), , drop = FALSE]
  }

  list(data_matrix = data_matrix, metadata = metadata)
}


#' Create a Volcano Plot with Optional Classification Colouring
#'
#' Significant lipids (adj.P.Val < \code{pval_threshold} AND
#' |logFC| > \code{logfc_threshold}) are coloured; non-significant points
#' are grey. When \code{classification_data} is supplied, significant lipids
#' are coloured by the selected classification column.
#'
#' @param results Data frame of differential analysis results (must have
#'   columns \code{logFC} and \code{adj.P.Val}; row names = lipid names).
#' @param title Plot title.
#' @param logfc_threshold Absolute log-fold-change threshold.
#' @param pval_threshold Adjusted p-value threshold.
#' @param top_labels Number of top significant lipids to label.
#' @param classification_data Optional classification data frame with a
#'   \code{Lipid} column.
#' @param color_by Column in \code{classification_data} to use for colouring.
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' res <- perform_differential_analysis(norm, d$metadata, "Sample Group",
#'   contrasts_list = NULL, method = "limma"
#' )
#' p <- create_volcano_plot_labeled(res$results[[1]])
#' print(p)
create_volcano_plot_labeled <- function(results,
                                        title = "Volcano Plot",
                                        logfc_threshold = 1,
                                        pval_threshold = 0.05,
                                        top_labels = 15,
                                        classification_data = NULL,
                                        color_by = NULL) {
  if (!is.data.frame(results)) results <- as.data.frame(results)
  if (!all(c("logFC", "adj.P.Val") %in% colnames(results))) {
    stop(
      "`results` must contain \"logFC\" and \"adj.P.Val\" columns ",
      "(e.g. one contrast's data frame from perform_differential_analysis()$results)."
    )
  }
  if (!"Lipid" %in% colnames(results)) results$Lipid <- rownames(results)

  results$Significance <- "Not Significant"
  results$Significance[results$adj.P.Val < pval_threshold &
    results$logFC > logfc_threshold] <- "Increased"
  results$Significance[results$adj.P.Val < pval_threshold &
    results$logFC < -logfc_threshold] <- "Decreased"
  results$IsSignificant <- results$Significance != "Not Significant"
  results$neg_log10_pval <- -log10(results$adj.P.Val)

  results <- results[order(results$adj.P.Val), ]
  top_features <- utils::head(results$Lipid[results$IsSignificant], top_labels)
  results$Label <- ifelse(results$Lipid %in% top_features, results$Lipid, "")

  if (!is.null(classification_data) && !is.null(color_by)) {
    # Deduplicate classification by Lipid (keep first row per lipid name)
    # to avoid inflating the results data frame after the merge.
    classification_data <- classification_data[
      !duplicated(classification_data$Lipid), ,
      drop = FALSE
    ]

    results <- merge(results, classification_data, by = "Lipid", all.x = TRUE)

    # Check whether the requested column exists after the merge
    no_match_msg <- NULL
    if (!color_by %in% colnames(results)) {
      no_match_msg <- paste0(
        "Column '", color_by,
        "' not found after merge. Check classification file columns."
      )
    } else {
      n_matched <- sum(!is.na(results[[color_by]]) & results$IsSignificant)
      if (n_matched == 0 && any(results$IsSignificant)) {
        ex_cls <- if (nrow(classification_data) > 0) {
          classification_data$Lipid[1]
        } else {
          "none"
        }
        ex_data <- results$Lipid[results$IsSignificant][1]
        no_match_msg <- paste0(
          "No lipid names matched between classification and data for column '",
          color_by, "'. ",
          "Classification example: '", ex_cls, "' | ",
          "Data example: '", ex_data, "'. ",
          "Lipid names must match exactly (case-sensitive)."
        )
      }
    }

    if (!is.null(no_match_msg)) {
      # Fall back to significance coloring and add subtitle warning
      p <- ggplot2::ggplot(
        results,
        ggplot2::aes(
          x = logFC, y = neg_log10_pval,
          color = Significance
        )
      ) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::scale_color_manual(
          values = c(
            Increased = "red", Decreased = "blue",
            "Not Significant" = "grey50"
          ),
          name = "Abundance Change"
        ) +
        ggplot2::labs(subtitle = paste("WARNING:", no_match_msg))
    } else {
      results$ColorGroup <- ifelse(
        results$IsSignificant,
        as.character(results[[color_by]]),
        "Not Significant"
      )
      results$ColorGroup[is.na(results$ColorGroup)] <- "Unknown"

      unique_classes <- unique(results$ColorGroup[results$ColorGroup != "Not Significant"])
      unique_classes <- unique_classes[!is.na(unique_classes)]
      class_colors <- if (length(unique_classes) > 0) {
        stats::setNames(scales::hue_pal()(length(unique_classes)), unique_classes)
      } else {
        c()
      }
      all_colors <- c(class_colors, "Not Significant" = "grey70")

      p <- ggplot2::ggplot(
        results,
        ggplot2::aes(
          x = logFC, y = neg_log10_pval,
          color = ColorGroup
        )
      ) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::scale_color_manual(
          values = all_colors,
          name   = color_by,
          breaks = c(unique_classes, "Not Significant")
        )
    }
  } else {
    p <- ggplot2::ggplot(
      results,
      ggplot2::aes(
        x = logFC, y = neg_log10_pval,
        color = Significance
      )
    ) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::scale_color_manual(
        values = c(
          Increased = "red", Decreased = "blue",
          "Not Significant" = "grey50"
        ),
        name = "Abundance Change"
      )
  }

  p +
    ggplot2::geom_vline(
      xintercept = c(-logfc_threshold, logfc_threshold),
      linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pval_threshold),
      linetype = "dashed", color = "grey40"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Label),
      max.overlaps = 20, size = 3, box.padding = 0.5, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}


#' Create a Lipid Expression Barplot Ordered by Group
#'
#' Produces per-lipid barplots coloured by sample group. Samples are
#' automatically sorted by group (then alphabetically within group) for a
#' cleaner visual.
#'
#' @param data_matrix Numeric matrix (samples as rows, lipids as columns).
#' @param metadata Metadata data frame.
#' @param selected_lipids Character vector of lipid names to plot.
#' @param selected_samples Optional character vector of sample names to retain.
#' @param selected_groups Optional character vector of group names to retain.
#' @param group_column Name of the group column in \code{metadata}.
#' @param data_type Label for the y-axis subtitle ("raw" or "normalized").
#' @return A single \code{ggplot2} object (one lipid) or a named list of
#'   \code{ggplot2} objects (multiple lipids).
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' p <- create_lipid_expression_barplot(
#'   d$numeric_data, d$metadata,
#'   selected_lipids = colnames(d$numeric_data)[1],
#'   group_column = "Sample Group"
#' )
#' print(p)
create_lipid_expression_barplot <- function(data_matrix,
                                            metadata,
                                            selected_lipids,
                                            selected_samples = NULL,
                                            selected_groups = NULL,
                                            group_column = "Sample Group",
                                            data_type = "normalized") {
  .validate_metadata(metadata)
  if (missing(selected_lipids) || length(selected_lipids) == 0) {
    stop("`selected_lipids` must be a non-empty character vector of lipid names.")
  }
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)
  # Orient with samples as rows, resolved by sample identity rather than by
  # comparing dimension sizes (which fails when n_samples > n_features).
  data_matrix <- .orient_matrix(data_matrix, metadata, want = "samples_rows")
  if (is.null(rownames(data_matrix))) {
    rownames(data_matrix) <- paste0("Sample_", seq_len(nrow(data_matrix)))
  }

  metadata_sample_names <- if ("Sample Name" %in% colnames(metadata)) {
    metadata$`Sample Name`
  } else {
    rn <- rownames(metadata)
    if (is.null(rn)) paste0("Sample_", seq_len(nrow(metadata))) else rn
  }

  # Filter by groups
  if (!is.null(selected_groups) && length(selected_groups) > 0 &&
    group_column %in% colnames(metadata)) {
    group_mask <- metadata[[group_column]] %in% selected_groups
    metadata <- metadata[group_mask, , drop = FALSE]
    metadata_sample_names <- metadata_sample_names[group_mask]
  }

  # Filter by individual samples
  if (!is.null(selected_samples) && length(selected_samples) > 0) {
    sample_mask <- metadata_sample_names %in% selected_samples
    metadata <- metadata[sample_mask, , drop = FALSE]
    metadata_sample_names <- metadata_sample_names[sample_mask]
  }

  # Align data matrix to metadata
  common_samples <- intersect(rownames(data_matrix), metadata_sample_names)
  if (length(common_samples) > 0) {
    data_matrix <- data_matrix[common_samples, , drop = FALSE]
    if ("Sample Name" %in% colnames(metadata)) {
      metadata <- metadata[match(common_samples, metadata$`Sample Name`), ,
        drop = FALSE
      ]
    } else {
      metadata <- metadata[match(common_samples, rownames(metadata)), ,
        drop = FALSE
      ]
    }
  } else {
    n_s <- min(nrow(data_matrix), nrow(metadata))
    data_matrix <- data_matrix[seq_len(n_s), , drop = FALSE]
    metadata <- metadata[seq_len(n_s), , drop = FALSE]
  }

  if (nrow(data_matrix) == 0 || nrow(metadata) == 0) {
    stop("No samples remaining after filtering. Check sample/group selection.")
  }

  available_lipids <- intersect(selected_lipids, colnames(data_matrix))
  if (length(available_lipids) == 0) {
    stop("None of the selected lipids are present in the data.")
  }

  data_subset <- data_matrix[, available_lipids, drop = FALSE]

  # Pre-compute a group-sorted sample order (used for all lipid plots)
  group_vec <- if (group_column %in% colnames(metadata)) {
    as.character(metadata[[group_column]])
  } else {
    rep("Unknown", nrow(metadata))
  }
  sample_order <- rownames(data_subset)[order(group_vec, rownames(data_subset))]

  plots <- lapply(available_lipids, function(lipid) {
    plot_data <- data.frame(
      Sample = rownames(data_subset),
      Abundance = as.numeric(data_subset[, lipid]),
      Group = group_vec,
      stringsAsFactors = FALSE
    )
    # Apply group-sorted factor for x-axis
    plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = Sample, y = Abundance, fill = Group)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(
        title    = paste("Expression of", lipid),
        subtitle = paste("Data type:", data_type),
        x        = "Sample",
        y        = "Abundance"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom"
      ) +
      ggplot2::scale_fill_brewer(palette = "Set2")
  })

  names(plots) <- available_lipids
  if (length(plots) == 1) plots[[1]] else plots
}


#' Create a PCA Plot with Optional Confidence or Visual Ellipses
#'
#' The colour/fill legends are merged so that ellipses do not introduce
#' duplicate legend keys. Sample labels are kept separate from group labels.
#'
#' @param pca_data Data frame with columns \code{PC1}, \code{PC2}, \code{Group}
#'   and (optionally) \code{Sample} for point labels.
#' @param variance_explained Numeric vector of length \eqn{\ge 2} with the
#'   percent variance explained by PC1 and PC2.
#' @param ellipse_type One of \code{"none"}, \code{"confidence"}, or
#'   \code{"visual"}.
#' @param confidence_level Numeric confidence level for \code{"confidence"}
#'   ellipses (default \code{0.95}).
#' @param title Plot title.
#' @param show_sample_labels Logical. If \code{TRUE}, sample names are shown
#'   as text labels next to each point.
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' pca_res <- perform_pca(norm, d$metadata, "Sample Group")
#' p <- create_pca_plot_with_ellipses(pca_res$pca_data, pca_res$variance_explained)
#' print(p)
create_pca_plot_with_ellipses <- function(pca_data,
                                          variance_explained,
                                          ellipse_type = "none",
                                          confidence_level = 0.95,
                                          title = "PCA Analysis",
                                          show_sample_labels = FALSE) {
  if (!is.data.frame(pca_data) || !all(c("PC1", "PC2", "Group") %in% colnames(pca_data))) {
    stop(
      "`pca_data` must be a data.frame with \"PC1\", \"PC2\", and \"Group\" ",
      "columns (e.g. the `$pca_data` element of perform_pca())."
    )
  }
  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = Group)
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.7) +
    ggplot2::labs(
      title = title,
      x     = paste0("PC1 (", round(variance_explained[1], 2), "% variance)"),
      y     = paste0("PC2 (", round(variance_explained[2], 2), "% variance)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Add sample name labels using geom_text (plotly-compatible).
  # ggrepel::geom_text_repel is intentionally NOT used here because plotly
  # does not support it and produces a warning. geom_text with check_overlap
  # gives a clean result in both static and plotly contexts.
  if (show_sample_labels && "Sample" %in% colnames(pca_data)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = Sample),
      size = 2.5,
      vjust = -0.8,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  # Add ellipses and merge colour+fill legends to avoid duplicate keys
  if (ellipse_type == "confidence") {
    p <- p +
      ggplot2::stat_ellipse(
        ggplot2::aes(fill = Group),
        type = "norm",
        level = confidence_level,
        geom = "polygon",
        alpha = 0.1,
        color = NA # border drawn by the line layer below
      ) +
      ggplot2::stat_ellipse(
        type = "norm",
        level = confidence_level,
        linetype = 2,
        linewidth = 1,
        show.legend = FALSE
      ) +
      # Merge colour and fill into one legend key
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(fill = NA)),
        fill  = "none"
      )
  } else if (ellipse_type == "visual") {
    p <- p +
      ggplot2::stat_ellipse(
        ggplot2::aes(fill = Group),
        type = "t",
        level = 0.8,
        geom = "polygon",
        alpha = 0.1,
        color = NA
      ) +
      ggplot2::stat_ellipse(
        type = "t",
        level = 0.8,
        linetype = 2,
        linewidth = 1,
        show.legend = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(fill = NA)),
        fill  = "none"
      )
  }

  p
}


#' Create a PLS-DA Plot with Optional Ellipses
#'
#' @param plsda_data Data frame with columns \code{Comp1}, \code{Comp2},
#'   \code{Group}, and (optionally) \code{Sample}.
#' @param ellipse_type One of \code{"none"}, \code{"confidence"}, or
#'   \code{"visual"}.
#' @param confidence_level Numeric confidence level (default \code{0.95}).
#' @param title Plot title.
#' @param show_sample_labels Logical. Show sample name labels if \code{TRUE}.
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' res <- perform_plsda(norm, d$metadata, "Sample Group")
#' p <- create_plsda_plot_with_ellipses(res$scores_data)
#' print(p)
create_plsda_plot_with_ellipses <- function(plsda_data,
                                            ellipse_type = "none",
                                            confidence_level = 0.95,
                                            title = "PLS-DA Analysis",
                                            show_sample_labels = FALSE) {
  if (!is.data.frame(plsda_data) || !all(c("Comp1", "Comp2", "Group") %in% colnames(plsda_data))) {
    stop(
      "`plsda_data` must be a data.frame with \"Comp1\", \"Comp2\", and \"Group\" ",
      "columns (e.g. the `$scores_data` element of perform_plsda())."
    )
  }
  p <- ggplot2::ggplot(
    plsda_data,
    ggplot2::aes(x = Comp1, y = Comp2, color = Group)
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::labs(
      title = title,
      x = "PLS-DA Component 1",
      y = "PLS-DA Component 2"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  if (show_sample_labels && "Sample" %in% colnames(plsda_data)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = Sample),
      size = 2.5, vjust = -0.8,
      check_overlap = TRUE, show.legend = FALSE
    )
  }

  if (ellipse_type == "confidence") {
    p <- p +
      ggplot2::stat_ellipse(ggplot2::aes(fill = Group),
        type = "norm", level = confidence_level,
        geom = "polygon", alpha = 0.1, color = NA
      ) +
      ggplot2::stat_ellipse(
        type = "norm", level = confidence_level,
        linetype = 2, linewidth = 1, show.legend = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(fill = NA)),
        fill = "none"
      )
  } else if (ellipse_type == "visual") {
    p <- p +
      ggplot2::stat_ellipse(ggplot2::aes(fill = Group),
        type = "t", level = 0.8,
        geom = "polygon", alpha = 0.1, color = NA
      ) +
      ggplot2::stat_ellipse(
        type = "t", level = 0.8,
        linetype = 2, linewidth = 1, show.legend = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(fill = NA)),
        fill = "none"
      )
  }

  p
}


#' Create an Enrichment Dotplot
#'
#' @param enrichment_data Data frame of fgsea results.
#' @param title Plot title.
#' @param max_pathways Maximum number of pathways (top by p-value).
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' cls <- classify_lipids(colnames(norm))
#' res <- perform_differential_analysis(norm, d$metadata, "Sample Group",
#'   contrasts_list = NULL, method = "limma"
#' )
#' enrich <- perform_enrichment_analysis(res$results, cls, min_set_size = 3)
#' p <- create_enrichment_dotplot(enrich[[1]][["LipidGroup"]])
#' print(p)
create_enrichment_dotplot <- function(enrichment_data,
                                      title = "Enrichment Analysis",
                                      max_pathways = 15) {
  if (!is.data.frame(enrichment_data)) {
    stop("`enrichment_data` must be a data.frame, not ", class(enrichment_data)[1], ".")
  }
  if (nrow(enrichment_data) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = "No enrichment results"
      ) +
      ggplot2::theme_void())
  }

  enrichment_data <- utils::head(
    enrichment_data[order(enrichment_data$pval), ],
    max_pathways
  )

  ggplot2::ggplot(
    enrichment_data,
    ggplot2::aes(
      x = NES,
      y = stats::reorder(pathway, NES),
      size = size,
      color = -log10(padj)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(
      low = "blue", high = "red",
      name = "-log10(padj)"
    ) +
    ggplot2::scale_size_continuous(name = "Set Size") +
    ggplot2::labs(
      title = title,
      x = "Normalized Enrichment Score", y = "Pathway"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )
}


#' Create an Enrichment Barplot
#'
#' @param enrichment_data Data frame of fgsea results.
#' @param title Plot title.
#' @param max_pathways Maximum number of pathways (top by p-value).
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' cls <- classify_lipids(colnames(norm))
#' res <- perform_differential_analysis(norm, d$metadata, "Sample Group",
#'   contrasts_list = NULL, method = "limma"
#' )
#' enrich <- perform_enrichment_analysis(res$results, cls, min_set_size = 3)
#' p <- create_enrichment_barplot(enrich[[1]][["LipidGroup"]])
#' print(p)
create_enrichment_barplot <- function(enrichment_data,
                                      title = "Enrichment Analysis",
                                      max_pathways = 15) {
  if (!is.data.frame(enrichment_data)) {
    stop("`enrichment_data` must be a data.frame, not ", class(enrichment_data)[1], ".")
  }
  if (nrow(enrichment_data) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text",
        x = 0.5, y = 0.5,
        label = "No enrichment results"
      ) +
      ggplot2::theme_void())
  }

  enrichment_data <- utils::head(
    enrichment_data[order(enrichment_data$pval), ],
    max_pathways
  )
  enrichment_data$Direction <- ifelse(
    enrichment_data$NES > 0, "Increased abundance", "Decreased abundance"
  )

  ggplot2::ggplot(
    enrichment_data,
    ggplot2::aes(
      x = NES,
      y = stats::reorder(pathway, NES),
      fill = Direction
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(
      values = c("Increased abundance" = "red", "Decreased abundance" = "blue")
    ) +
    ggplot2::labs(
      title = title,
      x = "Normalized Enrichment Score", y = "Pathway"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y     = ggplot2::element_text(size = 8),
      legend.position = "bottom"
    )
}
