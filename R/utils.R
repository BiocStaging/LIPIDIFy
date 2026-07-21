#' Load Lipidomics Data
#'
#' Reads a CSV file and separates metadata columns from numeric lipid abundance columns.
#'
#' @param file_path Path to the lipidomics data file (CSV format).
#' @param metadata_columns Character vector of expected metadata column names.
#' @return A named list with components: \code{data} (original data frame),
#'   \code{metadata} (metadata data frame), \code{numeric_data} (numeric matrix).
#' @export
#' @examples
#' # Write a minimal CSV then load it
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(
#'     "Sample Name" = c("S1", "S2"), "Sample Group" = c("A", "B"),
#'     "PC 16:0" = c(1000, 2000), check.names = FALSE
#'   ),
#'   tmp,
#'   row.names = FALSE
#' )
#' loaded <- load_lipidomics_data(tmp)
#' dim(loaded$numeric_data)
#' unlink(tmp)
load_lipidomics_data <- function(file_path,
                                 metadata_columns = c(
                                   "Sample Name", "Sample Group",
                                   "Tumour ID", "Weight (mg)"
                                 )) {
  data <- utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

  # Remove QC pool samples if present
  if ("Sample Group" %in% colnames(data)) {
    data <- data[!grepl("PBQC", data$`Sample Group`, ignore.case = TRUE), ]
  }

  # Identify metadata vs. numeric columns
  metadata_cols <- names(data)[names(data) %in% metadata_columns]
  candidate_numeric_cols <- names(data)[!names(data) %in% metadata_columns]

  actually_numeric <- vapply(candidate_numeric_cols, function(col) {
    converted <- tryCatch(
      as.numeric(data[[col]]),
      warning = function(w) rep(NA_real_, length(data[[col]])),
      error = function(e) rep(NA_real_, length(data[[col]]))
    )
    !all(is.na(converted))
  }, logical(1))

  non_numeric_cols <- candidate_numeric_cols[!actually_numeric]
  numeric_cols <- candidate_numeric_cols[actually_numeric]
  all_metadata_cols <- unique(c(metadata_cols, non_numeric_cols))

  metadata <- data[, names(data) %in% all_metadata_cols, drop = FALSE]
  numeric_data <- data[, numeric_cols, drop = FALSE]

  for (col in names(numeric_data)) {
    numeric_data[[col]] <- as.numeric(numeric_data[[col]])
  }

  all_na_cols <- vapply(numeric_data, function(x) all(is.na(x)), logical(1))
  if (any(all_na_cols)) {
    numeric_data <- numeric_data[, !all_na_cols, drop = FALSE]
  }

  numeric_data <- as.matrix(numeric_data)
  storage.mode(numeric_data) <- "numeric"

  if ("Sample Name" %in% colnames(metadata)) {
    rownames(numeric_data) <- metadata$`Sample Name`
    rownames(metadata) <- metadata$`Sample Name`
  }

  # Drop a trailing column that is entirely NA or constant
  if (ncol(numeric_data) > 0) {
    last_col <- numeric_data[, ncol(numeric_data)]
    if (all(is.na(last_col)) ||
      length(unique(last_col[!is.na(last_col)])) <= 1) {
      numeric_data <- numeric_data[, -ncol(numeric_data), drop = FALSE]
    }
  }

  list(data = data, metadata = metadata, numeric_data = numeric_data)
}

# ---------------------------------------------------------------------------
# Matrix orientation
# ---------------------------------------------------------------------------

#' Orient a Data Matrix to Samples- or Features-as-Rows (internal)
#'
#' Deterministically orients an abundance matrix without guessing from the
#' relative number of samples and features. Orientation is resolved by matching
#' the matrix dimension names against the sample identifiers in \code{metadata}
#' (the \code{Sample Name} column, or the metadata row names). This is robust
#' whether there are more samples than features or vice versa. When sample
#' identifiers cannot be matched by name, it falls back to matching the sample
#' \emph{count}; if that is also ambiguous it assumes the matrix already has
#' samples as rows.
#'
#' @param data_matrix Numeric matrix or object coercible to one.
#' @param metadata Optional metadata data frame used to identify samples.
#' @param want Either \code{"samples_rows"} or \code{"features_rows"}.
#' @return \code{data_matrix} transposed as needed so that the requested entity
#'   is on the rows.
#' @keywords internal
.orient_matrix <- function(data_matrix, metadata = NULL,
                           want = c("samples_rows", "features_rows")) {
  want <- match.arg(want)
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)

  sample_ids <- NULL
  n_samples <- NA_integer_
  if (!is.null(metadata) && nrow(metadata) > 0) {
    n_samples <- nrow(metadata)
    if ("Sample Name" %in% colnames(metadata)) {
      sample_ids <- as.character(metadata[["Sample Name"]])
    } else if (!is.null(rownames(metadata))) {
      sample_ids <- rownames(metadata)
    }
  }

  rows_are_samples <- NA
  # 1) Preferred: match sample identifiers against dimension names
  if (!is.null(sample_ids)) {
    if (!is.null(rownames(data_matrix)) &&
      length(intersect(rownames(data_matrix), sample_ids)) > 0) {
      rows_are_samples <- TRUE
    } else if (!is.null(colnames(data_matrix)) &&
      length(intersect(colnames(data_matrix), sample_ids)) > 0) {
      rows_are_samples <- FALSE
    }
  }
  # 2) Fallback: match the sample count (only when unambiguous)
  if (is.na(rows_are_samples) && !is.na(n_samples)) {
    if (nrow(data_matrix) == n_samples && ncol(data_matrix) != n_samples) {
      rows_are_samples <- TRUE
    } else if (ncol(data_matrix) == n_samples && nrow(data_matrix) != n_samples) {
      rows_are_samples <- FALSE
    }
  }
  # 3) Last resort: assume samples are already on the rows
  if (is.na(rows_are_samples)) rows_are_samples <- TRUE

  if (want == "samples_rows" && !rows_are_samples) {
    data_matrix <- t(data_matrix)
  } else if (want == "features_rows" && rows_are_samples) {
    data_matrix <- t(data_matrix)
  }
  data_matrix
}

# ---------------------------------------------------------------------------
# Lipid classification
# ---------------------------------------------------------------------------

#' Classify Lipids Based on Their Names
#'
#' Classifies a vector of lipid names into lipid group, type, and saturation
#' category using regular-expression pattern matching.
#'
#' @param lipid_names Character vector of lipid names.
#' @return Data frame with columns \code{Lipid}, \code{LipidGroup},
#'   \code{LipidType}, and \code{Saturation}.
#' @export
#' @examples
#' classify_lipids(c("PC 16:0_18:1", "TG 16:0_18:1_20:4", "Cer 16:0"))
classify_lipids <- function(lipid_names) {
  classify_single <- function(lipid) {
    lc <- stringr::str_trim(lipid)

    group <- "Other/Unclassified"
    type <- "Other"
    saturation <- "Unclassified"

    classes <- list(
      list(pat = "^PC", grp = "Glycerophospholipids", typ = "Phosphatidylcholine"),
      list(pat = "^PE", grp = "Glycerophospholipids", typ = "Phosphatidylethanolamine"),
      list(pat = "^PG", grp = "Glycerophospholipids", typ = "Phosphatidylglycerol"),
      list(pat = "^PI", grp = "Glycerophospholipids", typ = "Phosphatidylinositol"),
      list(pat = "^PS", grp = "Glycerophospholipids", typ = "Phosphatidylserine"),
      list(pat = "^LPC", grp = "Lysophospholipids", typ = "Lysophosphatidylcholine"),
      list(pat = "^LPE", grp = "Lysophospholipids", typ = "Lysophosphatidylethanolamine"),
      list(pat = "^LPI", grp = "Lysophospholipids", typ = "Lysophosphatidylinositol"),
      list(pat = "^DG", grp = "Glycerolipids", typ = "Diacylglycerol"),
      list(pat = "^TG", grp = "Glycerolipids", typ = "Triacylglycerol"),
      list(pat = "^dhCer", grp = "Sphingolipids", typ = "Dihydroceramide"),
      list(pat = "^Cer", grp = "Sphingolipids", typ = "Ceramide"),
      list(pat = "^Hex", grp = "Sphingolipids", typ = "Hexosylceramide"),
      list(pat = "^SM", grp = "Sphingolipids", typ = "Sphingomyelin"),
      list(pat = "^GM", grp = "Sphingolipids", typ = "Ganglioside"),
      list(pat = "^Sulfatide", grp = "Sphingolipids", typ = "Sulfatide"),
      list(pat = "^CE", grp = "Sterol Lipids", typ = "Cholesteryl Ester"),
      list(pat = "^Desmosterol", grp = "Sterol Lipids", typ = "Desmosterol"),
      list(pat = "^COH", grp = "Sterol Lipids", typ = "Cholesterol/Related"),
      list(pat = "^Ubiquinone", grp = "Sterol Lipids", typ = "Ubiquinone"),
      list(pat = "^AcylCarnitine", grp = "Acylcarnitines", typ = "Acylcarnitine"),
      list(pat = "(IS)", grp = "Internal Standard", typ = lc)
    )

    for (cl in classes) {
      if (stringr::str_detect(lc, stringr::regex(cl$pat, ignore_case = TRUE))) {
        group <- cl$grp
        type <- cl$typ
        break
      }
    }

    saturation <- determine_saturation(lc)

    data.frame(
      Lipid = lc, LipidGroup = group, LipidType = type,
      Saturation = saturation, stringsAsFactors = FALSE
    )
  }

  do.call(rbind, lapply(lipid_names, classify_single))
}


#' Determine Fatty-Acid Saturation from a Lipid Name
#'
#' Parses standard lipid name notation to count double bonds and classifies the
#' species as SFA (0 double bonds), MUFA (1) or PUFA (>1).
#'
#' @param lipid_name A single character string with the lipid name.
#' @return One of \code{"SFA"}, \code{"MUFA"}, \code{"PUFA"}, or
#'   \code{"Unclassified"}.
#' @export
#' @examples
#' determine_saturation("PC 16:0_18:1") # "MUFA"
#' determine_saturation("PE 18:0_18:0") # "SFA"
#' determine_saturation("TG 16:0_18:1_20:4") # "PUFA"
determine_saturation <- function(lipid_name) {
  classify_db <- function(db) {
    if (db == 0) {
      return("SFA")
    }
    if (db == 1) {
      return("MUFA")
    }
    if (db > 1) {
      return("PUFA")
    }
    "Unclassified"
  }

  tryCatch(
    {
      # Pattern 1 - standard "16:0" tokens
      # Note: \b word boundaries are intentionally omitted here.
      # Underscore is a word character in regex, so \b would NOT match
      # the boundary before "18" in "16:0_18:1", causing only the first
      # fatty acid to be captured and giving a wrong SFA result.
      p1 <- stringr::str_extract_all(lipid_name, "\\d+:\\d+")[[1]]
      if (length(p1) > 0) {
        db <- sum(vapply(p1, function(x) as.numeric(strsplit(x, ":")[[1]][2]), numeric(1)),
          na.rm = TRUE
        )
        return(classify_db(db))
      }

      # Pattern 2 - parenthesis notation e.g. "PE(P-18:0_20:5)"
      p2 <- stringr::str_extract_all(lipid_name, "\\(.*?\\)")[[1]]
      if (length(p2) > 0) {
        inner <- stringr::str_remove_all(p2[1], "[()]")
        fa <- stringr::str_extract_all(inner, "\\d+:\\d+")[[1]]
        if (length(fa) > 0) {
          db <- sum(vapply(fa, function(x) as.numeric(strsplit(x, ":")[[1]][2]), numeric(1)),
            na.rm = TRUE
          )
          return(classify_db(db))
        }
      }

      # Pattern 3 - space-separated " 18:2"
      p3 <- stringr::str_extract_all(lipid_name, "\\s\\d+:\\d+")[[1]]
      if (length(p3) > 0) {
        db <- sum(
          vapply(
            stringr::str_trim(p3),
            function(x) as.numeric(strsplit(x, ":")[[1]][2]), numeric(1)
          ),
          na.rm = TRUE
        )
        return(classify_db(db))
      }

      # Pattern 4 - total notation "PE 35:2" (two-/three-digit carbon count at end)
      p4 <- stringr::str_extract(lipid_name, "\\s(\\d{2,3}):(\\d+)\\s*$")
      if (!is.na(p4)) {
        m <- stringr::str_match(p4, "\\s\\d{2,3}:(\\d+)")
        if (!is.na(m[1, 2])) {
          return(classify_db(as.numeric(m[1, 2])))
        }
      }

      "Unclassified"
    },
    error = function(e) "Unclassified"
  )
}


#' Test Saturation Classification
#'
#' Convenience function to verify saturation detection on a set of lipid names.
#'
#' @param test_lipids Optional character vector of lipid names to test.
#'   If \code{NULL}, a built-in test set is used.
#' @return A data frame with columns \code{Lipid} and \code{Saturation},
#'   printed to the console and returned invisibly.
#' @export
#' @examples
#' test_saturation_classification()
test_saturation_classification <- function(test_lipids = NULL) {
  if (is.null(test_lipids)) {
    test_lipids <- c(
      "PC 16:0_18:1", # MUFA  (0+1=1)
      "LPC 18:2", # PUFA  (2)
      "PE 18:0_18:0", # SFA   (0+0=0)
      "TG 16:0_18:1_20:4", # PUFA  (0+1+4=5)
      "PE(P-18:0_20:5)", # PUFA  (0+5=5)
      "SM 16:0", # SFA   (0)
      "CE 18:1", # MUFA  (1)
      "PE 35:2", # PUFA  (total 2)
      "PC 34:1", # MUFA  (total 1)
      "LPC 20:5" # PUFA  (5)
    )
  }

  results <- data.frame(
    Lipid = test_lipids,
    Saturation = vapply(test_lipids, determine_saturation, character(1)),
    stringsAsFactors = FALSE
  )

  message("Saturation Classification Test Results:")
  message(paste(utils::capture.output(results), collapse = "\n"))
  invisible(results)
}

# ---------------------------------------------------------------------------
# Normalization methods
# ---------------------------------------------------------------------------

#' Return Available Normalization Method Names
#'
#' @return Character vector of normalization method names supported by
#'   \code{apply_normalizations}.
#' @export
#' @examples
#' get_normalization_methods()
get_normalization_methods <- function() {
  c(
    "TIC", "PQN", "Quantile", "Log2Median", "Median", "Mean",
    "Log2", "Log10", "Sqrt", "None"
  )
}

#' Return Human-Readable Descriptions of Normalization Methods
#'
#' Used by the Shiny app to populate help text.
#'
#' @return Named character vector (name = method key, value = description).
#' @export
#' @examples
#' descs <- get_normalization_descriptions()
#' cat(descs["TIC"])
get_normalization_descriptions <- function() {
  c(
    TIC = paste(
      "Total Ion Current (TIC) normalization divides each sample by its total",
      "signal intensity, then rescales to the global mean TIC.",
      "Best suited when the total amount of material injected varies between samples."
    ),
    PQN = paste(
      "Probabilistic Quotient Normalization (PQN) computes a reference spectrum",
      "(median across samples), calculates per-sample quotients, and normalizes",
      "by the median quotient. More robust to large fold-changes than TIC."
    ),
    Quantile = paste(
      "Quantile normalization forces all samples to have an identical",
      "intensity distribution by replacing each value with the mean of",
      "identically ranked values across samples.",
      "After this step, per-sample boxplots will look nearly identical",
      "-- this is expected and correct behaviour, not a bug."
    ),
    Log2Median = paste(
      "Log2 Median Centering: log2-transforms the data (log2(x + 1)),",
      "then subtracts each sample's median deviation from the global median,",
      "effectively centering samples on a common baseline.",
      "Useful when variance scales with the mean intensity.",
      "Note: this is a simplified variance-stabilising step, not the full",
      "maximum-likelihood VSN procedure (Huber et al. 2002)."
    ),
    Median = paste(
      "Median normalization scales each sample so that its median equals the",
      "global median across all samples. Simple and robust to outliers."
    ),
    Mean = paste(
      "Mean normalization scales each sample so that its mean equals the global",
      "mean. Similar to Median normalization but sensitive to extreme values."
    ),
    Log2 = "Log base-2 transformation. Standard for mass-spectrometry data before statistical analysis.",
    Log10 = "Log base-10 transformation. Reduces right-skew; useful for data spanning several orders of magnitude.",
    Sqrt = "Square-root transformation. A milder variance-stabilising alternative to log transforms.",
    None = "No normalization applied. Data passed through unchanged."
  )
}

#' Apply a Sequence of Normalization Methods
#'
#' Applies normalization methods in the order supplied, passing the output of
#' each step as the input to the next.
#'
#' @param data Numeric matrix with samples in rows and lipids in columns.
#' @param methods Character vector of method names as returned by
#'   \code{get_normalization_methods()}.
#' @return Normalized numeric matrix of the same dimensions as \code{data}.
#' @export
#' @examples
#' m <- matrix(rlnorm(60, 8, 1), nrow = 6, ncol = 10)
#' apply_normalizations(m, c("TIC", "Log2"))
apply_normalizations <- function(data, methods) {
  if (!is.matrix(data)) data <- as.matrix(data)
  original_rownames <- rownames(data)

  normalized <- data
  for (method in methods) {
    normalized <- switch(method,
      TIC        = normalize_tic(normalized),
      PQN        = normalize_pqn(normalized),
      Quantile   = normalize_quantile(normalized),
      Log2Median = normalize_log2median(normalized),
      Median     = normalize_median(normalized),
      Mean     = normalize_mean(normalized),
      Log2     = log2(normalized + 1),
      Log10    = log10(normalized + 1),
      Sqrt     = sqrt(abs(normalized)),
      None     = normalized,
      normalized # unknown method: pass through
    )
  }

  if (!is.null(original_rownames)) rownames(normalized) <- original_rownames
  normalized
}

# ----- Individual normalization functions ----------------------------------

#' TIC Normalization
#'
#' Divides each sample by its total ion current (row sum) and rescales to the
#' global mean TIC.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return TIC-normalized matrix of the same dimensions.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_tic(m)
normalize_tic <- function(data) {
  tic <- rowSums(data, na.rm = TRUE)
  mean_tic <- mean(tic, na.rm = TRUE)
  sweep(data, 1, tic, "/") * mean_tic
}

#' PQN Normalization
#'
#' Probabilistic Quotient Normalization. Uses the per-feature median across
#' samples as the reference spectrum.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return PQN-normalized matrix.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_pqn(m)
normalize_pqn <- function(data) {
  reference <- apply(data, 2, stats::median, na.rm = TRUE)
  quotients <- sweep(data, 2, reference, "/")
  median_quotients <- apply(quotients, 1, stats::median, na.rm = TRUE)
  sweep(data, 1, median_quotients, "/")
}

#' Quantile Normalization
#'
#' Forces all samples to share an identical intensity distribution.
#' After this step per-sample boxplots will look nearly identical --
#' that is the intended and correct behaviour of quantile normalization.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return Quantile-normalized matrix of the same dimensions.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_quantile(m)
normalize_quantile <- function(data) {
  if (!is.matrix(data)) data <- as.matrix(data)
  storage.mode(data) <- "numeric"

  all_na_cols <- apply(data, 2, function(x) all(is.na(x)))
  if (any(all_na_cols)) data <- data[, !all_na_cols, drop = FALSE]

  if (nrow(data) < 2 || ncol(data) < 2) {
    warning(
      "Quantile normalization requires at least 2 samples and 2 features. ",
      "Returning input unchanged."
    )
    return(data)
  }

  # Transpose: features as rows, samples as columns
  data_t <- t(data)
  ranks <- apply(data_t, 2, rank, ties.method = "average")
  mean_sorted <- rowMeans(apply(data_t, 2, sort), na.rm = TRUE)
  normalized_t <- data_t
  for (i in seq_len(ncol(data_t))) {
    normalized_t[, i] <- mean_sorted[ranks[, i]]
  }
  t(normalized_t)
}

#' Log2 Median Centering Normalization
#'
#' Log2-transforms the data (log2(x + 1)), then median-centres each sample
#' relative to the global median. This is a simplified variance-stabilising
#' step suitable for mass-spectrometry lipidomics data.
#'
#' @details
#' This method is \strong{not} equivalent to the full VSN procedure of
#' Huber et al. (2002), which uses maximum-likelihood estimation. If true
#' VSN is required, use the \pkg{vsn} Bioconductor package directly.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return Log2-median-centred numeric matrix.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_log2median(m)
normalize_log2median <- function(data) {
  log_data <- log2(data + 1)
  medians <- apply(log_data, 1, stats::median, na.rm = TRUE)
  global_median <- stats::median(medians, na.rm = TRUE)
  normalized <- sweep(log_data, 1, medians - global_median, "-")
  2^normalized - 1
}

#' @rdname normalize_log2median
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_vsn(m)  # deprecated alias
normalize_vsn <- function(data) {
  .Deprecated("normalize_log2median",
    msg = paste0(
      "'normalize_vsn' has been renamed to 'normalize_log2median' to better ",
      "reflect its implementation (log2 + median centering). ",
      "Please update your code."
    )
  )
  normalize_log2median(data)
}

#' Median Normalization
#'
#' Scales each sample so its median equals the global median across all samples.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return Median-normalized matrix.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_median(m)
normalize_median <- function(data) {
  medians <- apply(data, 1, stats::median, na.rm = TRUE)
  global_median <- stats::median(medians, na.rm = TRUE)
  sweep(data, 1, medians / global_median, "/")
}

#' Mean Normalization
#'
#' Scales each sample so its mean equals the global mean across all samples.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @return Mean-normalized matrix.
#' @export
#' @examples
#' m <- matrix(c(1000, 2000, 3000, 4000, 500, 1500), nrow = 2)
#' normalize_mean(m)
normalize_mean <- function(data) {
  means <- apply(data, 1, mean, na.rm = TRUE)
  global_mean <- mean(means, na.rm = TRUE)
  sweep(data, 1, means / global_mean, "/")
}


#' Normalize Lipidomics Data
#'
#' Convenience wrapper around \code{apply_normalizations}.
#'
#' @param data Numeric matrix (samples as rows, lipids as columns).
#' @param methods Character vector of normalization method names.
#' @return Normalized numeric matrix.
#' @export
#' @examples
#' m <- matrix(rlnorm(60, 8, 1), nrow = 6, ncol = 10)
#' normalize_lipidomics_data(m, c("TIC", "Log2"))
normalize_lipidomics_data <- function(data, methods = c("TIC", "Log2")) {
  apply_normalizations(data, methods)
}

# ---------------------------------------------------------------------------
# Custom classification helpers
# ---------------------------------------------------------------------------

#' Load Custom Lipid Classification from a CSV File
#'
#' @param file_path Path to a CSV file with a \code{Lipid} column followed by
#'   one or more classification columns.
#' @return Data frame with \code{Lipid} as the first column.
#' @export
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(
#'     Lipid = c("PC 16:0_18:1", "PE 18:0"),
#'     Class = c("Phospholipid", "Phospholipid")
#'   ),
#'   tmp,
#'   row.names = FALSE
#' )
#' cls <- load_custom_classification(tmp)
#' unlink(tmp)
load_custom_classification <- function(file_path) {
  classification <- utils::read.csv(file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"Lipid" %in% colnames(classification)) {
    stop("Classification file must contain a 'Lipid' column.")
  }
  other_cols <- classification[, colnames(classification) != "Lipid", drop = FALSE]
  data.frame(Lipid = classification$Lipid, other_cols, stringsAsFactors = FALSE)
}

#' Export Lipid Classification to CSV
#'
#' @param classification Data frame with lipid classifications.
#' @param file_path Output file path.
#' @return Invisibly returns \code{TRUE} on success.
#' @export
#' @examples
#' cls <- classify_lipids(c("PC 16:0_18:1", "PE 18:0_20:4"))
#' tmp <- tempfile(fileext = ".csv")
#' export_classification(cls, tmp)
#' unlink(tmp)
export_classification <- function(classification, file_path) {
  utils::write.csv(classification, file_path, row.names = FALSE)
  invisible(TRUE)
}

#' Get Lipid Classification
#'
#' Wrapper around \code{classify_lipids} for convenient scripting use.
#'
#' @param lipid_names Character vector of lipid names.
#' @return Data frame with columns \code{Lipid}, \code{LipidGroup},
#'   \code{LipidType}, and \code{Saturation}.
#' @export
#' @examples
#' get_lipid_classification(c("PC 16:0_18:1", "PE 18:0_20:4"))
get_lipid_classification <- function(lipid_names) {
  classify_lipids(lipid_names)
}

#' Load Custom Enrichment Sets from a CSV File
#'
#' The CSV must have columns \code{Lipid} and \code{Set_Name}.
#' A lipid may appear in multiple rows to belong to multiple sets.
#'
#' @param file_path Path to the CSV file.
#' @return Named list of character vectors (one vector per set).
#' @export
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(
#'   data.frame(
#'     Lipid = c("PC 16:0_18:1", "PE 18:0"),
#'     Set_Name = c("Phospholipids", "Phospholipids")
#'   ),
#'   tmp,
#'   row.names = FALSE
#' )
#' sets <- load_custom_enrichment_sets(tmp)
#' unlink(tmp)
load_custom_enrichment_sets <- function(file_path) {
  sets_df <- utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("Lipid", "Set_Name") %in% colnames(sets_df))) {
    stop("Custom sets file must contain 'Lipid' and 'Set_Name' columns.")
  }
  split(sets_df$Lipid, sets_df$Set_Name)
}

# ---------------------------------------------------------------------------
# Missing value imputation
# ---------------------------------------------------------------------------

#' Return Available Imputation Method Names
#'
#' @return Character vector of imputation method names supported by
#'   \code{impute_missing_values}.
#' @export
#' @examples
#' get_imputation_methods()
get_imputation_methods <- function() {
  c("half_min", "min", "zero", "mean", "median", "knn")
}

#' Return Human-Readable Descriptions of Imputation Methods
#'
#' @return Named character vector (name = method key, value = description).
#' @export
#' @examples
#' descs <- get_imputation_descriptions()
#' cat(descs["half_min"])
get_imputation_descriptions <- function() {
  c(
    half_min = paste(
      "Half-minimum: replaces each NA with half the column minimum.",
      "Models the limit of detection (LOD) and is the most widely used",
      "approach in untargeted MS lipidomics."
    ),
    min = paste(
      "Minimum: replaces each NA with the column minimum observed value.",
      "Slightly more conservative than half-minimum."
    ),
    zero = paste(
      "Zero: replaces NA with 0.",
      "Use only when absence of signal genuinely means zero abundance."
    ),
    mean = paste(
      "Mean: replaces each NA with the column mean of observed values.",
      "Simple but can introduce bias when data are not missing at random (MNAR)."
    ),
    median = paste(
      "Median: replaces each NA with the column median of observed values.",
      "More robust than mean to extreme values."
    ),
    knn = paste(
      "K-Nearest Neighbours (KNN): imputes each missing value using the",
      "k most similar samples (Euclidean distance in feature space).",
      "Requires the impute Bioconductor package.",
      "Recommended when data are missing at random (MAR)."
    )
  )
}

#' Impute Missing Values in a Lipidomics Data Matrix
#'
#' Replaces \code{NA} values using the chosen strategy. Imputation should be
#' applied \strong{before} normalisation so that missing values do not bias
#' per-sample scaling factors.
#'
#' @param data_matrix Numeric matrix with samples as rows and lipids as
#'   columns.
#' @param method Imputation method. One of \code{"half_min"} (default),
#'   \code{"min"}, \code{"zero"}, \code{"mean"}, \code{"median"},
#'   \code{"knn"}. See \code{\link{get_imputation_descriptions}} for details.
#' @param k Integer. Number of nearest neighbours for \code{"knn"}
#'   (ignored for other methods). Default \code{5}.
#' @param seed Integer. Random seed for reproducibility when \code{method =
#'   "knn"}. Default \code{42}.
#' @return Imputed numeric matrix of the same dimensions as
#'   \code{data_matrix}.
#' @export
#' @examples
#' m <- matrix(c(1000, NA, 3000, NA, 500, 1500), nrow = 2)
#' impute_missing_values(m, method = "half_min")
impute_missing_values <- function(data_matrix,
                                   method = "half_min",
                                   k      = 5L,
                                   seed   = 42L) {
  method <- match.arg(method, get_imputation_methods())
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"

  n_missing <- sum(is.na(data_matrix))
  if (n_missing == 0L) {
    message("No missing values detected. Returning data unchanged.")
    return(data_matrix)
  }
  message(sprintf(
    "Imputing %d missing value(s) (%.1f%% of matrix) using method: %s",
    n_missing,
    100 * n_missing / length(data_matrix),
    method
  ))

  result <- switch(method,
    half_min = {
      apply(data_matrix, 2L, function(col) {
        col[is.na(col)] <- min(col, na.rm = TRUE) / 2
        col
      })
    },
    min = {
      apply(data_matrix, 2L, function(col) {
        col[is.na(col)] <- min(col, na.rm = TRUE)
        col
      })
    },
    zero = {
      data_matrix[is.na(data_matrix)] <- 0
      data_matrix
    },
    mean = {
      apply(data_matrix, 2L, function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
    },
    median = {
      apply(data_matrix, 2L, function(col) {
        col[is.na(col)] <- stats::median(col, na.rm = TRUE)
        col
      })
    },
    knn = {
      if (!requireNamespace("impute", quietly = TRUE)) {
        warning(
          "Package 'impute' is required for KNN imputation. ",
          "Install with: BiocManager::install('impute'). ",
          "Falling back to half-minimum imputation.",
          call. = FALSE
        )
        return(impute_missing_values(data_matrix, method = "half_min"))
      }
      set.seed(seed)
      # impute::impute.knn expects features as rows, samples as columns
      imp <- impute::impute.knn(t(data_matrix), k = k)
      t(imp$data)
    }
  )

  rownames(result) <- rownames(data_matrix)
  colnames(result) <- colnames(data_matrix)
  storage.mode(result) <- "numeric"
  result
}


# ---------------------------------------------------------------------------
# Batch effect correction
# ---------------------------------------------------------------------------

#' Correct Batch Effects from a Normalised Lipidomics Matrix
#'
#' Removes known technical batch effects while preserving biological signal.
#' Batch correction should be applied \strong{after} normalisation.
#'
#' Two methods are supported:
#' \describe{
#'   \item{\code{"limma"}}{Uses \code{\link[limma]{removeBatchEffect}}.
#'     Requires only the limma package (already a dependency). Suitable for
#'     most experimental designs.}
#'   \item{\code{"combat"}}{Uses \code{sva::ComBat} with parametric
#'     empirical Bayes adjustment. Requires the \pkg{sva} Bioconductor
#'     package (\code{BiocManager::install("sva")}). Generally more robust
#'     when batch effects are large.}
#' }
#'
#' @param data_matrix Numeric matrix with samples as rows and lipids as
#'   columns (typically the output of \code{\link{apply_normalizations}}).
#' @param metadata Data frame of sample metadata aligned with
#'   \code{data_matrix} (same row order).
#' @param batch_column Character. Name of the column in \code{metadata}
#'   containing batch labels.
#' @param group_column Character. Name of the biological group column to
#'   protect from removal (default \code{"Sample Group"}). Pass \code{NULL}
#'   to skip group protection.
#' @param method One of \code{"limma"} (default) or \code{"combat"}.
#' @return Batch-corrected numeric matrix of the same dimensions as
#'   \code{data_matrix}.
#' @export
#' @examples
#' d    <- load_lipidomics_data_from_df(generate_example_data())
#' norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
#' d$metadata$Batch <- rep(c("Batch1", "Batch2"), each = 10)
#' corrected <- correct_batch_effects(norm, d$metadata,
#'                                    batch_column = "Batch")
#' dim(corrected)
correct_batch_effects <- function(data_matrix,
                                   metadata,
                                   batch_column,
                                   group_column = "Sample Group",
                                   method       = "limma") {
  method <- match.arg(method, c("limma", "combat"))
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)

  if (!batch_column %in% colnames(metadata)) {
    stop(
      "batch_column '", batch_column, "' not found in metadata. ",
      "Available columns: ", paste(colnames(metadata), collapse = ", "),
      call. = FALSE
    )
  }

  batch <- as.factor(as.character(metadata[[batch_column]]))

  # Build a design matrix to protect biological signal.
  # If the group structure is perfectly confounded with batch (i.e. the design
  # matrix is rank-deficient after including batch), limma::removeBatchEffect
  # will throw a cryptic error. We detect this and fall back to correcting
  # without group protection, emitting a warning instead.
  design_protect <- NULL
  if (!is.null(group_column) && group_column %in% colnames(metadata)) {
    groups         <- factor(make.names(as.character(metadata[[group_column]])))
    design_raw     <- stats::model.matrix(~ groups)

    # Check rank: combine batch indicator columns with the group design and
    # test whether the full matrix is full-rank.
    batch_dummy  <- stats::model.matrix(~ batch - 1)
    combined     <- cbind(design_raw, batch_dummy)
    if (qr(combined)$rank < ncol(combined)) {
      warning(
        "The group column ('", group_column, "') is perfectly confounded with ",
        "the batch column ('", batch_column, "'). ",
        "Proceeding without group protection -- all between-group variance may ",
        "be removed. Consider a design where batches contain samples from ",
        "multiple groups.",
        call. = FALSE
      )
      design_protect <- NULL
    } else {
      design_protect <- design_raw
    }
  }

  # limma and ComBat expect features x samples. Resolve orientation by sample
  # identity first (robust to either input layout), then transpose so that
  # features are on the rows and samples on the columns.
  data_matrix <- .orient_matrix(data_matrix, metadata, want = "samples_rows")
  data_t <- t(data_matrix)

  corrected_t <- switch(method,
    limma = {
      limma::removeBatchEffect(data_t,
        batch  = batch,
        design = design_protect
      )
    },
    combat = {
      if (!requireNamespace("sva", quietly = TRUE)) {
        warning(
          "Package 'sva' is required for ComBat. ",
          "Install with: BiocManager::install('sva'). ",
          "Falling back to limma::removeBatchEffect.",
          call. = FALSE
        )
        return(correct_batch_effects(data_matrix, metadata, batch_column,
                                      group_column, method = "limma"))
      }
      sva::ComBat(
        dat        = data_t,
        batch      = batch,
        mod        = design_protect,
        par.prior  = TRUE,
        prior.plots = FALSE
      )
    }
  )

  corrected <- t(corrected_t)
  rownames(corrected) <- rownames(data_matrix)
  colnames(corrected) <- colnames(data_matrix)
  corrected
}

