#' Generate Example Dataset
#'
#' Creates a synthetic lipidomics dataset with realistic lipid names and group differences
#' suitable for demonstrating differential analysis capabilities.
#'
#' @return Data frame with simulated lipidomics data including 4 groups with 5 replicates each
#' @export
#' @examples
#' example_data <- generate_example_data()
#' head(example_data[, 1:10])
generate_example_data <- function() {
  # Sample metadata with generic group names
  sample_names <- paste0("Sample_", seq_len(20))
  sample_groups <- rep(c("GroupA", "GroupB", "GroupC", "GroupD"), each = 5)
  tumor_ids <- rep(seq_len(10), each = 2)
  weights <- stats::rnorm(20, mean = 50, sd = 10)

  # Lipid features (simulate realistic lipid names and intensities)
  lipid_classes <- c("PC", "PE", "PI", "PS", "PG", "LPC", "LPE", "LPI", "SM", "Cer", "TG", "DG", "CE")
  n_lipids_per_class <- c(50, 40, 20, 3, 4, 30, 10, 4, 25, 25, 50, 15, 25)

  lipid_names <- c()
  for (i in seq_along(lipid_classes)) {
    class_name <- lipid_classes[i]
    n_lipids <- n_lipids_per_class[i]

    # Generate realistic lipid names
    for (j in seq_len(n_lipids)) {
      # Random fatty acid combinations
      fa1 <- sample(c(14:22), 1)
      fa2 <- sample(c(14:24), 1)
      db1 <- sample(c(0:2), 1)
      db2 <- sample(c(0:6), 1)

      if (class_name %in% c("LPC", "LPE", "LPI")) {
        lipid_name <- paste0(class_name, " ", fa1, ":", db1)
      } else if (class_name %in% c("Cer", "SM")) {
        lipid_name <- paste0(class_name, " ", fa1, ":", db1)
      } else {
        lipid_name <- paste0(class_name, " ", fa1, ":", db1, "_", fa2, ":", db2)
      }

      lipid_names <- c(lipid_names, lipid_name)
    }
  }

  # Make lipid names unique
  lipid_names <- make.unique(lipid_names)

  # Generate intensity data with stronger group effects for better differential analysis
  n_features <- length(lipid_names)
  intensity_data <- matrix(nrow = 20, ncol = n_features)
  colnames(intensity_data) <- lipid_names
  rownames(intensity_data) <- sample_names

  # Base intensities
  for (i in seq_len(n_features)) {
    base_intensity <- stats::rlnorm(20, meanlog = 8, sdlog = 1)

    # Add stronger group effects for more features
    if (i %% 8 == 0) { # GroupB upregulated
      group_effect <- ifelse(sample_groups == "GroupB", 3.0, 1)
      intensity_data[, i] <- base_intensity * group_effect
    } else if (i %% 10 == 0) { # GroupC upregulated
      group_effect <- ifelse(sample_groups == "GroupC", 2.5, 1)
      intensity_data[, i] <- base_intensity * group_effect
    } else if (i %% 12 == 0) { # GroupD downregulated
      group_effect <- ifelse(sample_groups == "GroupD", 0.4, 1)
      intensity_data[, i] <- base_intensity * group_effect
    } else if (i %% 15 == 0) { # GroupA vs others
      group_effect <- ifelse(sample_groups == "GroupA", 2.0, 1)
      intensity_data[, i] <- base_intensity * group_effect
    } else {
      intensity_data[, i] <- base_intensity
    }
  }

  # Combine metadata and intensity data
  example_data <- data.frame(
    "Sample Name" = sample_names,
    "Sample Group" = sample_groups,
    "Tumour ID" = tumor_ids,
    "Weight (mg)" = weights,
    intensity_data,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  return(example_data)
}

#' Load Lipidomics Data from Data Frame
#'
#' Processes a lipidomics data frame by separating metadata and numeric data columns.
#'
#' @param data_df Data frame with lipidomics data
#' @param metadata_columns Vector of metadata column names
#' @return List containing data components (data, metadata, numeric_data)
#' @export
#' @examples
#' data_df <- generate_example_data()
#' loaded_data <- load_lipidomics_data_from_df(data_df)
#' names(loaded_data)
load_lipidomics_data_from_df <- function(data_df, metadata_columns = c("Sample Name", "Sample Group", "Tumour ID", "Weight (mg)")) {
  if (!is.data.frame(data_df)) {
    stop("`data_df` must be a data.frame, not ", class(data_df)[1], ".")
  }

  # Filter out PBQC samples if present
  if ("Sample Group" %in% colnames(data_df)) {
    data_df <- data_df[!grepl("PBQC", data_df$`Sample Group`, ignore.case = TRUE), ]
  }

  # Start with user-specified metadata columns
  metadata_cols <- names(data_df)[names(data_df) %in% metadata_columns]

  # Identify candidate numeric columns (those NOT in metadata_columns)
  candidate_numeric_cols <- names(data_df)[!names(data_df) %in% metadata_columns]

  # Also move any non-numeric columns into metadata (fallback by column type)
  actually_numeric <- vapply(candidate_numeric_cols, function(col) {
    converted <- tryCatch(
      as.numeric(data_df[[col]]),
      warning = function(w) rep(NA_real_, length(data_df[[col]])),
      error = function(e) rep(NA_real_, length(data_df[[col]]))
    )
    !all(is.na(converted))
  }, logical(1))

  # Columns that fail numeric coercion go to metadata
  non_numeric_cols <- candidate_numeric_cols[!actually_numeric]
  numeric_cols <- candidate_numeric_cols[actually_numeric]

  # Combine all metadata columns
  all_metadata_cols <- unique(c(metadata_cols, non_numeric_cols))

  # Separate metadata and numeric data
  metadata <- data_df[, names(data_df) %in% all_metadata_cols, drop = FALSE]
  numeric_data <- data_df[, numeric_cols, drop = FALSE]

  # Explicitly coerce each numeric column to numeric
  for (col in names(numeric_data)) {
    numeric_data[[col]] <- as.numeric(numeric_data[[col]])
  }

  # Drop columns that are all NA after coercion
  all_na_cols <- vapply(
    seq_len(ncol(numeric_data)),
    function(j) all(is.na(numeric_data[, j])), logical(1)
  )
  if (any(all_na_cols)) {
    numeric_data <- numeric_data[, !all_na_cols, drop = FALSE]
  }

  # Convert to matrix only after ensuring all columns are numeric
  numeric_data <- as.matrix(numeric_data)
  storage.mode(numeric_data) <- "numeric"

  # Set rownames based on Sample Name
  if ("Sample Name" %in% colnames(metadata)) {
    rownames(numeric_data) <- metadata$`Sample Name`
    rownames(metadata) <- metadata$`Sample Name`
  }

  # Remove the last column if it appears to be irrelevant
  if (ncol(numeric_data) > 0) {
    last_col <- numeric_data[, ncol(numeric_data)]
    if (all(is.na(last_col)) || length(unique(last_col[!is.na(last_col)])) <= 1) {
      numeric_data <- numeric_data[, -ncol(numeric_data), drop = FALSE]
    }
  }

  return(list(
    data = data_df,
    metadata = metadata,
    numeric_data = numeric_data
  ))
}
