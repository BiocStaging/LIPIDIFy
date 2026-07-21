#' Perform Differential Analysis
#'
#' @param data_matrix Normalized data matrix (features as rows, samples as columns)
#' @param metadata Metadata data frame. Must contain \code{group_column} and
#'   have a number of rows matching the samples in \code{data_matrix}.
#' @param group_column Column name in metadata containing group information
#' @param contrasts_list List of contrasts to perform
#' @param method Method to use: "limma" or "edger"
#' @param design Optional user-supplied design matrix (samples as rows). If
#'   \code{NULL} (the default), a design matrix of the form \code{~0 + group}
#'   is built automatically from \code{group_column}. When supplying a custom
#'   design, \code{contrasts_list} must also be supplied (limma/edgeR-style
#'   contrast strings referencing the design's column names).
#' @return List containing results for each contrast
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' res <- perform_differential_analysis(norm, d$metadata, "Sample Group",
#'   contrasts_list = NULL, method = "limma"
#' )
#' names(res)
perform_differential_analysis <- function(data_matrix, metadata, group_column = "Sample Group",
                                          contrasts_list = NULL, method = "limma",
                                          design = NULL) {
  if (!method %in% c("limma", "edger")) {
    stop("`method` must be either \"limma\" or \"edger\", not \"", method, "\".")
  }
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data.frame, not ", class(metadata)[1], ".")
  }
  if (!group_column %in% colnames(metadata)) {
    stop(
      "`group_column` (\"", group_column, "\") was not found in `metadata`. ",
      "Available columns: ", paste(colnames(metadata), collapse = ", ")
    )
  }
  if (!is.null(design) && is.null(contrasts_list)) {
    stop("`contrasts_list` must be supplied when a custom `design` matrix is used.")
  }

  if (method == "limma") {
    .perform_differential_analysis_limma(data_matrix, metadata, group_column, contrasts_list, design)
  } else {
    .perform_differential_analysis_edger(data_matrix, metadata, group_column, contrasts_list, design)
  }
}

#' Perform Differential Analysis with limma
#'
#' @param data_matrix Normalized data matrix (features as rows, samples as columns)
#' @param metadata Metadata data frame
#' @param group_column Column name in metadata containing group information
#' @param contrasts_list List of contrasts to perform
#' @param design Optional user-supplied design matrix; see
#'   \code{\link{perform_differential_analysis}}.
#' @return List containing LIMMA results for each contrast
#' @noRd
.perform_differential_analysis_limma <- function(data_matrix, metadata, group_column = "Sample Group",
                                                 contrasts_list = NULL, design = NULL) {
  # Convert to matrix if needed
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }

  # Ensure data matrix is oriented with features as rows (samples as columns),
  # resolved by sample identity rather than by comparing dimension sizes.
  data_matrix <- .orient_matrix(data_matrix, metadata, want = "features_rows")

  # Store original group labels for UI/display, and a make.names()-sanitized
  # version for use as limma design/contrast identifiers.
  original_groups <- as.character(metadata[[group_column]])
  original_levels <- unique(original_groups)

  sanitized_groups <- make.names(original_groups)
  sanitized_levels <- make.names(original_levels)

  level_mapping <- sanitized_levels
  names(level_mapping) <- original_levels
  reverse_mapping <- original_levels
  names(reverse_mapping) <- sanitized_levels

  if (is.null(design)) {
    groups <- factor(sanitized_groups, levels = sanitized_levels)
    design <- stats::model.matrix(~ 0 + groups)
    colnames(design) <- sanitized_levels
  } else {
    design <- .validate_design(design, data_matrix)
    # A user-supplied design is used as-is; contrast strings are expected to
    # reference its column names directly, so skip level sanitisation/mapping.
    level_mapping <- NULL
    reverse_mapping <- NULL
  }

  # Fit linear model
  fit <- limma::lmFit(data_matrix, design)

  # If no contrasts specified, create default contrasts (using sanitized names)
  if (is.null(contrasts_list)) {
    contrasts_list <- create_default_contrasts(sanitized_levels)
  } else if (!is.null(level_mapping)) {
    # Convert user-provided contrasts (written in terms of the original group
    # labels) to the sanitized names used as design column names.
    contrasts_list <- vapply(contrasts_list, function(contrast) {
      for (orig in original_levels) {
        sanitized <- level_mapping[[orig]]
        # Replace whole word matches only
        contrast <- gsub(
          paste0("\\b", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", orig), "\\b"),
          sanitized, contrast
        )
      }
      return(contrast)
    }, character(1), USE.NAMES = FALSE)
  }

  # Create contrast matrix
  contrast_matrix <- limma::makeContrasts(contrasts = contrasts_list, levels = design)

  # Fit contrasts and apply empirical Bayes
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  # Extract results for each contrast
  results <- list()
  original_contrast_names <- colnames(contrast_matrix)

  for (i in seq_len(ncol(contrast_matrix))) {
    # Convert sanitized contrast names back to the original group labels for
    # display, when a default (auto-sanitized) design was used.
    contrast_name <- original_contrast_names[i]
    display_name <- contrast_name
    if (!is.null(reverse_mapping)) {
      for (sanitized in sanitized_levels) {
        orig <- reverse_mapping[[sanitized]]
        display_name <- gsub(paste0("\\b", sanitized, "\\b"), orig, display_name)
      }
    }

    results[[display_name]] <- limma::topTable(fit2, coef = i, number = Inf, adjust.method = "fdr")
  }

  return(list(
    fit = fit2,
    results = results,
    design = design,
    contrasts = contrast_matrix,
    method = "limma",
    level_mapping = level_mapping,
    reverse_mapping = reverse_mapping
  ))
}

#' Perform Differential Analysis with EdgeR
#'
#' @param data_matrix Normalized data matrix (features as rows, samples as columns)
#' @param metadata Metadata data frame
#' @param group_column Column name in metadata containing group information
#' @param contrasts_list List of contrasts to perform
#' @param design Optional user-supplied design matrix; see
#'   \code{\link{perform_differential_analysis}}.
#' @return List containing EdgeR results for each contrast
#' @noRd
.perform_differential_analysis_edger <- function(data_matrix, metadata, group_column = "Sample Group",
                                                 contrasts_list = NULL, design = NULL) {
  # Statistical note: edgeR's negative-binomial model was developed for integer
  # count data (RNA-seq). Applying it to continuous mass-spectrometry lipidomics
  # intensities is an approximation. limma is the recommended method for
  # MS-based lipidomics. edgeR is provided here as an exploratory alternative;
  # results should be interpreted with caution and validated against limma.
  warning(
    "edgeR is designed for integer count data (e.g., RNA-seq). ",
    "Applying it to continuous lipidomics intensities is an approximation. ",
    "limma is the recommended method for MS-based lipidomics.",
    call. = FALSE
  )

  # Convert to matrix if needed
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }

  # Ensure data matrix is oriented with features as rows (samples as columns),
  # resolved by sample identity rather than by comparing dimension sizes.
  data_matrix <- .orient_matrix(data_matrix, metadata, want = "features_rows")

  # Shift any negative values to zero (can arise after log-scale normalisation)
  # without rounding to integers, preserving the continuous nature of the data.
  data_matrix[data_matrix < 0] <- 0

  # Store original group labels for UI/display, and a make.names()-sanitized
  # version for use as edgeR design/contrast identifiers.
  original_groups <- as.character(metadata[[group_column]])
  original_levels <- unique(original_groups)

  sanitized_groups <- make.names(original_groups)
  sanitized_levels <- make.names(original_levels)

  level_mapping <- sanitized_levels
  names(level_mapping) <- original_levels
  reverse_mapping <- original_levels
  names(reverse_mapping) <- sanitized_levels

  groups <- factor(sanitized_groups, levels = sanitized_levels)

  # Create DGEList object
  dge <- edgeR::DGEList(counts = data_matrix, group = groups)

  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge)

  if (is.null(design)) {
    design <- stats::model.matrix(~ 0 + groups)
    colnames(design) <- sanitized_levels
  } else {
    design <- .validate_design(design, data_matrix)
    # A user-supplied design is used as-is; contrast strings are expected to
    # reference its column names directly, so skip level sanitisation/mapping.
    level_mapping <- NULL
    reverse_mapping <- NULL
  }

  # If no contrasts specified, create default contrasts (using sanitized names)
  if (is.null(contrasts_list)) {
    contrasts_list <- create_default_contrasts(sanitized_levels)
  } else if (!is.null(level_mapping)) {
    # Convert user-provided contrasts (written in terms of the original group
    # labels) to the sanitized names used as design column names.
    contrasts_list <- vapply(contrasts_list, function(contrast) {
      for (orig in original_levels) {
        sanitized <- level_mapping[[orig]]
        # Replace whole word matches only
        contrast <- gsub(
          paste0("\\b", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", orig), "\\b"),
          sanitized, contrast
        )
      }
      return(contrast)
    }, character(1), USE.NAMES = FALSE)
  }

  # Create contrast matrix
  contrast_matrix <- limma::makeContrasts(contrasts = contrasts_list, levels = design)

  # Fit GLM
  fit <- edgeR::glmQLFit(dge, design)

  # Extract results for each contrast
  results <- list()
  original_contrast_names <- colnames(contrast_matrix)

  for (i in seq_len(ncol(contrast_matrix))) {
    # Convert sanitized contrast names back to the original group labels for
    # display, when a default (auto-sanitized) design was used.
    contrast_name <- original_contrast_names[i]
    display_name <- contrast_name
    if (!is.null(reverse_mapping)) {
      for (sanitized in sanitized_levels) {
        orig <- reverse_mapping[[sanitized]]
        display_name <- gsub(paste0("\\b", sanitized, "\\b"), orig, display_name)
      }
    }

    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix[, i])

    # Extract results in similar format to limma
    res <- edgeR::topTags(qlf, n = Inf)$table

    # Rename columns to match limma output
    colnames(res)[colnames(res) == "logCPM"] <- "AveExpr"
    colnames(res)[colnames(res) == "PValue"] <- "P.Value"
    colnames(res)[colnames(res) == "FDR"] <- "adj.P.Val"

    results[[display_name]] <- res
  }

  return(list(
    fit = fit,
    results = results,
    design = design,
    contrasts = contrast_matrix,
    method = "edger",
    level_mapping = level_mapping,
    reverse_mapping = reverse_mapping
  ))
}

#' Create Default Contrasts
#'
#' Creates default pairwise contrasts from group levels. If the levels contain
#' special characters that are invalid for limma/edgeR, they should already
#' be sanitized before calling this function.
#'
#' @param group_levels Character vector of group level names (e.g., c("Control","Treatment","Resistant"))
#' @return Character vector of limma-style contrast strings (e.g., "Treatment - Control")
#' @examples
#' create_default_contrasts(c("A", "B", "C"))
#' @export
create_default_contrasts <- function(group_levels) {
  group_levels <- as.character(group_levels)
  group_levels <- unique(group_levels[!is.na(group_levels)])
  if (length(group_levels) < 2) {
    return(character(0))
  }

  # Callers are expected to pass already-sanitized levels; make.names() here
  # is a safeguard for direct usage with raw group labels.
  safe_levels <- make.names(group_levels)

  contrasts <- c()
  for (i in seq_len(length(safe_levels) - 1)) {
    for (j in (i + 1):length(safe_levels)) {
      contrasts <- c(contrasts, paste0(safe_levels[j], " - ", safe_levels[i]))
    }
  }
  contrasts
}


#' Perform Enrichment Analysis
#'
#' @param results_list List of differential analysis results
#' @param classification_data Lipid classification data frame
#' @param min_set_size Minimum pathway set size
#' @param max_set_size Maximum pathway set size
#' @param custom_sets Optional named list of custom lipid sets
#' @return List containing GSEA results
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' cls <- classify_lipids(colnames(norm))
#' res <- perform_differential_analysis(norm, d$metadata, "Sample Group",
#'   contrasts_list = NULL, method = "limma"
#' )
#' enrich <- perform_enrichment_analysis(res$results, cls, min_set_size = 3)
#' names(enrich)
perform_enrichment_analysis <- function(results_list, classification_data,
                                        min_set_size = 5, max_set_size = 500,
                                        custom_sets = NULL) {
  if (!is.list(results_list) || is.data.frame(results_list) || is.null(names(results_list))) {
    stop("`results_list` must be a named list of differential analysis result ",
      "data frames (e.g. the `$results` element of `perform_differential_analysis()`).")
  }
  if (!is.data.frame(classification_data)) {
    stop("`classification_data` must be a data.frame, not ", class(classification_data)[1], ".")
  }
  if (!"Lipid" %in% colnames(classification_data)) {
    stop("`classification_data` must contain a \"Lipid\" column.")
  }

  enrichment_results <- list()

  for (contrast_name in names(results_list)) {
    message("Processing enrichment for: ", contrast_name)

    results <- results_list[[contrast_name]]

    # Create ranked list
    ranked_vector <- results$logFC
    names(ranked_vector) <- rownames(results)
    ranked_vector <- sort(ranked_vector, decreasing = TRUE)

    # Merge with classification
    merged_data <- merge(
      data.frame(Lipid = names(ranked_vector), logFC = ranked_vector, stringsAsFactors = FALSE),
      classification_data,
      by = "Lipid",
      all.x = TRUE
    )

    # Create pathway sets from classification columns
    pathway_sets_list <- list()

    # Get all classification columns (excluding Lipid and logFC)
    class_columns <- setdiff(colnames(classification_data), "Lipid")

    for (class_col in class_columns) {
      if (class_col %in% colnames(merged_data)) {
        sets <- .create_pathway_sets(merged_data, class_col)
        if (length(sets) > 0) {
          pathway_sets_list[[class_col]] <- sets
        }
      }
    }

    # Add custom sets if provided
    if (!is.null(custom_sets)) {
      pathway_sets_list[["custom"]] <- custom_sets
    }

    # Run GSEA for each set of pathways
    gsea_results <- list()
    for (set_name in names(pathway_sets_list)) {
      gsea_res <- .run_fgsea_safe(pathway_sets_list[[set_name]], ranked_vector, min_set_size, max_set_size)
      if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
        gsea_results[[set_name]] <- gsea_res
      }
    }

    enrichment_results[[contrast_name]] <- gsea_results
  }

  return(enrichment_results)
}

#' Create Pathway Sets
#'
#' @param merged_data Data frame with lipids and classifications
#' @param classification_column Column name for classification
#' @return Named list of pathway sets
#' @noRd
.create_pathway_sets <- function(merged_data, classification_column) {
  grouped_data <- split(merged_data$Lipid, merged_data[[classification_column]])

  # Remove any NULL or empty groups
  pathway_sets <- grouped_data[vapply(grouped_data, length, integer(1)) > 0L]

  return(pathway_sets)
}

#' Run FGSEA with Error Handling
#'
#' @param pathway_sets Named list of pathway sets
#' @param ranked_vector Named numeric vector of ranked statistics
#' @param min_size Minimum set size
#' @param max_size Maximum set size
#' @return FGSEA results data frame or NULL if error
#' @noRd
.run_fgsea_safe <- function(pathway_sets, ranked_vector, min_size, max_size) {
  if (length(pathway_sets) == 0) {
    return(data.frame())
  }

  tryCatch(
    {
      # Add small random noise to avoid ties
      ranked_vector <- ranked_vector + stats::rnorm(length(ranked_vector), sd = 1e-6)

      gsea_results <- fgsea::fgseaMultilevel(
        pathways = pathway_sets,
        stats = ranked_vector,
        minSize = min_size,
        maxSize = max_size
      )

      # Convert list columns to strings
      gsea_results <- .convert_list_columns_to_strings(gsea_results)

      return(gsea_results)
    },
    error = function(e) {
      # Emit a single, specific warning (avoid stacking redundant warnings)
      if (grepl("corrupt", e$message, ignore.case = TRUE)) {
        warning("fgsea database appears corrupt. Try reinstalling fgsea.",
          call. = FALSE
        )
      } else {
        warning("FGSEA error: ", conditionMessage(e), call. = FALSE)
      }
      return(NULL)
    }
  )
}

#' Convert List Columns to Strings
#'
#' @param df Data frame with potential list columns
#' @return Data frame with list columns converted to strings
#' @noRd
.convert_list_columns_to_strings <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- vapply(df[[col]], function(x) {
        if (is.null(x)) NA_character_ else paste(x, collapse = ",")
      }, character(1))
    }
  }
  return(df)
}

#' Perform PCA Analysis
#'
#' @param data_matrix Data matrix (samples as rows)
#' @param metadata Metadata data frame
#' @param group_column Group column name
#' @return List containing PCA results and plot
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' pca_res <- perform_pca(norm, d$metadata, "Sample Group")
#' names(pca_res)
perform_pca <- function(data_matrix, metadata, group_column = "Sample Group") {
  .validate_metadata(metadata, group_column)

  # Orient with samples as rows, resolved by sample identity rather than by
  # assuming features always outnumber samples (which fails when n_samples >
  # n_features, e.g. large cohorts with a focused lipid panel).
  data_matrix <- .orient_matrix(data_matrix, metadata, want = "samples_rows")

  # Remove any columns with zero variance
  var_cols <- apply(data_matrix, 2, stats::var, na.rm = TRUE) > 0
  data_matrix <- data_matrix[, var_cols, drop = FALSE]

  # Perform PCA
  pca_results <- FactoMineR::PCA(data_matrix, scale.unit = TRUE, graph = FALSE)

  # Create PCA plot data
  pca_data <- data.frame(
    Sample = rownames(data_matrix),
    PC1 = pca_results$ind$coord[, 1],
    PC2 = pca_results$ind$coord[, 2],
    Group = metadata[[group_column]],
    stringsAsFactors = FALSE
  )

  # Create plot
  variance_explained <- pca_results$eig[seq_len(2), 2]

  pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 4, alpha = 0.7) +
    ggplot2::labs(
      title = "PCA Analysis",
      x = paste0("PC1 (", round(variance_explained[1], 2), "% variance)"),
      y = paste0("PC2 (", round(variance_explained[2], 2), "% variance)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(list(
    pca_results = pca_results,
    pca_data = pca_data,
    plot = pca_plot,
    variance_explained = variance_explained
  ))
}

#' Perform PLS-DA Analysis
#'
#' @param data_matrix Data matrix (samples as rows)
#' @param metadata Metadata data frame
#' @param group_column Group column name
#' @param n_comp Number of components
#' @return List containing PLS-DA results and plot
#' @export
#' @examples
#' d <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' res <- perform_plsda(norm, d$metadata, "Sample Group")
#' names(res)
perform_plsda <- function(data_matrix, metadata, group_column = "Sample Group", n_comp = 2) {
  .validate_metadata(metadata, group_column)

  tryCatch(
    {
      if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
      }

      # Orient with samples as rows, resolved by sample identity rather than by
      # comparing dimension sizes (which fails when n_samples > n_features).
      data_matrix <- .orient_matrix(data_matrix, metadata, want = "samples_rows")

      # Remove any columns with zero variance
      var_cols <- apply(data_matrix, 2, stats::var, na.rm = TRUE) > 0
      data_matrix <- data_matrix[, var_cols, drop = FALSE]

      message("PLS-DA data matrix dimensions: ", paste(dim(data_matrix), collapse = " x "))
      message("PLS-DA metadata dimensions: ", paste(dim(metadata), collapse = " x "))

      # Prepare data
      groups <- factor(metadata[[group_column]])

      # Check for sample correspondence
      if (nrow(data_matrix) != length(groups)) {
        stop(sprintf(
          "Number of samples in data matrix (%d) does not match metadata (%d)",
          nrow(data_matrix), length(groups)
        ))
      }

      # Ensure we have at least 2 groups
      if (length(levels(groups)) < 2) {
        stop("PLS-DA requires at least 2 groups")
      }

      # Create dummy variables for groups
      group_dummy <- stats::model.matrix(~ groups - 1)

      # Handle column naming for stats::model.matrix() output
      colnames(group_dummy) <- paste0("Group_", levels(groups))

      message("Group dummy matrix dimensions: ", paste(dim(group_dummy), collapse = " x "))
      message("Unique groups: ", paste(levels(groups), collapse = ", "))

      # Perform PLS-DA using pls package
      plsda_results <- pls::plsr(group_dummy ~ data_matrix, ncomp = n_comp, validation = "LOO")

      # Extract scores
      scores_data <- data.frame(
        Sample = rownames(data_matrix),
        Comp1 = plsda_results$scores[, 1],
        Comp2 = if (n_comp > 1) plsda_results$scores[, 2] else rep(0, nrow(plsda_results$scores)),
        Group = groups,
        stringsAsFactors = FALSE
      )

      # Create plot
      plsda_plot <- ggplot2::ggplot(scores_data, ggplot2::aes(x = Comp1, y = Comp2, color = Group)) +
        ggplot2::geom_point(size = 4, alpha = 0.8) +
        ggplot2::labs(
          title = "PLS-DA Analysis",
          x = "PLS-DA Component 1",
          y = "PLS-DA Component 2"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")

      message("PLS-DA completed successfully")

      return(list(
        plsda_results = plsda_results,
        scores_data = scores_data,
        plot = plsda_plot
      ))
    },
    error = function(e) {
      # Single warning covers both notification and error propagation
      warning("PLS-DA failed: ", conditionMessage(e), call. = FALSE)

      # Fallback: create simple plot showing error
      error_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text",
          x = 0.5, y = 0.5,
          label = paste("PLS-DA failed:", e$message),
          size = 5
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "PLS-DA Analysis (Error)")

      return(list(
        plsda_results = NULL,
        scores_data = data.frame(),
        plot = error_plot
      ))
    }
  )
}
