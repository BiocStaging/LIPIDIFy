# =============================================================================
# tests/testthat/test-vsn-normalization.R
#
# Coverage for the VSN normalization method (normalize_vsn(), backed by
# vsn::justvsn()) and for the guarantee that VSN and Log2Median remain two
# distinct methods.
#
# Run with:
#   devtools::test()
# =============================================================================

# A fixed, fully deterministic fixture: values come from a closed-form
# expression rather than the RNG, so no seed is set and the caller's random
# stream is never touched. Positive, intensity-like, with a mean-dependent
# spread. vsn requires >= 42 features (its minDataPointsPerStratum default),
# so the default is 8 samples x 60 lipids.
vsn_matrix <- function(n_samples = 8, n_lipids = 60) {
  m <- outer(
    seq_len(n_samples), seq_len(n_lipids),
    function(i, j) exp(6 + 0.05 * j + 0.30 * sin(i + j) + 0.10 * i)
  )
  m <- matrix(m, nrow = n_samples, ncol = n_lipids)
  rownames(m) <- paste0("Sample_", seq_len(n_samples))
  colnames(m) <- paste0("Lipid_", seq_len(n_lipids))
  m
}

# 1. Method registry ---------------------------------------------------------
testthat::test_that("VSN appears in get_normalization_methods()", {
  testthat::expect_true("VSN" %in% get_normalization_methods())
})

testthat::test_that("Log2Median and VSN are separate, independently listed choices", {
  methods <- get_normalization_methods()
  testthat::expect_true("Log2Median" %in% methods)
  testthat::expect_true("VSN" %in% methods)
  testthat::expect_false(identical(
    which(methods == "VSN"), which(methods == "Log2Median")
  ))
  # Both must be described, and neither description may present them as
  # interchangeable.
  descs <- get_normalization_descriptions()
  testthat::expect_true(all(c("VSN", "Log2Median") %in% names(descs)))
  testthat::expect_match(descs[["VSN"]], "vsn::justvsn", fixed = TRUE)
  testthat::expect_match(descs[["Log2Median"]], "not the maximum-likelihood VSN")
})

testthat::test_that("Log2Median still performs log2-plus-median normalization", {
  m <- vsn_matrix()
  expected <- local({
    log_data <- log2(m + 1)
    medians <- apply(log_data, 1, stats::median)
    global_median <- stats::median(medians)
    2^sweep(log_data, 1, medians - global_median, "-") - 1
  })
  testthat::expect_equal(normalize_log2median(m), expected)
  # ...and is not an alias for the vsn implementation.
  testthat::skip_if_not_installed("vsn")
  testthat::expect_false(isTRUE(all.equal(
    normalize_log2median(m), normalize_vsn(m),
    tolerance = 1e-6
  )))
})

# 2. Output shape and values -------------------------------------------------
testthat::test_that("normalize_vsn preserves dimensions, dimnames and orientation", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  r <- normalize_vsn(m)
  testthat::expect_true(is.matrix(r))
  testthat::expect_type(r, "double")
  testthat::expect_equal(dim(r), dim(m))
  testthat::expect_equal(rownames(r), rownames(m))
  testthat::expect_equal(colnames(r), colnames(m))
  # Samples stay on the rows (the package normalization contract)
  testthat::expect_equal(nrow(r), 8L)
  testthat::expect_equal(ncol(r), 60L)
})

testthat::test_that("normalize_vsn returns finite numeric values for valid positive data", {
  testthat::skip_if_not_installed("vsn")
  r <- normalize_vsn(vsn_matrix())
  testthat::expect_true(all(is.finite(r)))
})

testthat::test_that("normalize_vsn agrees with a direct vsn::justvsn() call", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  # justvsn works on features x samples; normalize_vsn transposes explicitly.
  direct <- t(vsn::justvsn(t(m)))
  testthat::expect_equal(normalize_vsn(m), direct)
})

testthat::test_that("normalize_vsn accepts a coercible data frame", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  df <- as.data.frame(m)
  testthat::expect_equal(normalize_vsn(df), normalize_vsn(m))
})

# 3. Dispatch through apply_normalizations ----------------------------------
testthat::test_that("apply_normalizations('VSN') dispatches to normalize_vsn", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  testthat::expect_equal(apply_normalizations(m, "VSN"), normalize_vsn(m))
})

testthat::test_that("apply_normalizations dispatches VSN and Log2Median differently", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  testthat::expect_false(isTRUE(all.equal(
    apply_normalizations(m, "VSN"),
    apply_normalizations(m, "Log2Median"),
    tolerance = 1e-6
  )))
})

# 4. Missing dependency ------------------------------------------------------
testthat::test_that("missing vsn produces a clear error and no fallback", {
  m <- vsn_matrix()
  testthat::local_mocked_bindings(
    .vsn_installed = function() FALSE,
    .package = "LIPIDIFy"
  )
  testthat::expect_error(
    normalize_vsn(m),
    "requires the Bioconductor package 'vsn'",
    fixed = TRUE
  )
  testthat::expect_error(
    normalize_vsn(m),
    "BiocManager::install('vsn')",
    fixed = TRUE
  )
  # apply_normalizations must propagate the error, not silently return
  # Log2Median (or any other) output.
  testthat::expect_error(
    apply_normalizations(m, "VSN"),
    "requires the Bioconductor package 'vsn'",
    fixed = TRUE
  )
  testthat::expect_error(apply_normalizations(m, c("TIC", "VSN")))
})

# 5. Input validation --------------------------------------------------------
testthat::test_that("non-numeric input produces a descriptive error", {
  testthat::skip_if_not_installed("vsn")
  chr <- matrix(as.character(seq_len(8 * 60)), nrow = 8, ncol = 60)
  testthat::expect_error(normalize_vsn(chr), "must contain only numeric values")
  testthat::expect_error(normalize_vsn("not a matrix"), "must be a matrix or data.frame")
  df <- as.data.frame(vsn_matrix())
  df$Lipid_1 <- as.character(df$Lipid_1)
  testthat::expect_error(normalize_vsn(df), "must contain only numeric values")
})

testthat::test_that("infinite input produces a descriptive error", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  m[2, 3] <- Inf
  testthat::expect_error(normalize_vsn(m), "infinite value")
  m[2, 3] <- -Inf
  testthat::expect_error(normalize_vsn(m), "infinite value")
})

testthat::test_that("all-missing samples or features produce a descriptive error", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  m_row <- m
  m_row[3, ] <- NA_real_
  testthat::expect_error(
    suppressMessages(normalize_vsn(m_row)),
    "sample\\(s\\) contain no finite values"
  )
  m_col <- m
  m_col[, 5] <- NA_real_
  testthat::expect_error(
    suppressMessages(normalize_vsn(m_col)),
    "feature\\(s\\)[[:space:]]*contain no finite values"
  )
})

testthat::test_that("minDataPointsPerStratum is honoured and passed to vsn::justvsn", {
  testthat::skip_if_not_installed("vsn")
  small <- vsn_matrix(n_lipids = 30)
  # Rejected by default (30 < 42)...
  testthat::expect_error(normalize_vsn(small), "at least 42 lipid features")
  # ...and accepted, with the argument forwarded, when lowered explicitly.
  r <- normalize_vsn(small, minDataPointsPerStratum = 20L)
  testthat::expect_equal(dim(r), dim(small))
  testthat::expect_equal(
    r,
    t(vsn::justvsn(t(small), minDataPointsPerStratum = 20L))
  )
})

testthat::test_that("a sample with too few finite observations is rejected", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  m[5, -1] <- NA_real_ # sample 5 keeps a single observed value
  testthat::expect_error(
    suppressMessages(normalize_vsn(m)),
    "at least 2 finite observations per sample"
  )
})

testthat::test_that("too few samples or features produce a descriptive error", {
  testthat::skip_if_not_installed("vsn")
  testthat::expect_error(
    normalize_vsn(vsn_matrix(n_samples = 1)),
    "at least 2 samples"
  )
  testthat::expect_error(
    normalize_vsn(vsn_matrix(n_lipids = 20)),
    "at least 42 lipid features"
  )
})

# 6. Explicit missing-value behaviour ---------------------------------------
testthat::test_that("NA values are reported, preserved in place and never imputed", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  m[1, 1] <- NA_real_
  m[4, 7] <- NA_real_

  # The user is told explicitly that missing values are kept, not imputed.
  testthat::expect_message(normalize_vsn(m), "not imputed")

  r <- suppressMessages(normalize_vsn(m))
  testthat::expect_equal(dim(r), dim(m))
  testthat::expect_true(is.na(r[1, 1]))
  testthat::expect_true(is.na(r[4, 7]))
  testthat::expect_equal(sum(is.na(r)), 2L)
  testthat::expect_true(all(is.finite(r[!is.na(r)])))
})

testthat::test_that("negative values warn instead of being silently altered", {
  testthat::skip_if_not_installed("vsn")
  m <- vsn_matrix()
  m[2, 2] <- -50
  testthat::expect_warning(normalize_vsn(m), "negative value")
  r <- suppressWarnings(normalize_vsn(m))
  # The input matrix itself is untouched by the call.
  testthat::expect_equal(m[2, 2], -50)
  testthat::expect_equal(dim(r), dim(m))
})

# 7. Shiny dispatch ----------------------------------------------------------
testthat::test_that("the Shiny normalization dropdown offers VSN and Log2Median separately", {
  ui_src <- paste(deparse(LIPIDIFy:::.ui_tab_normalization), collapse = "\n")
  # Both checkbox groups are populated from get_normalization_methods()
  testthat::expect_match(ui_src, "get_normalization_methods", fixed = TRUE)
  choices <- get_normalization_methods()
  testthat::expect_true(all(c("VSN", "Log2Median") %in% choices))
})

# Wraps the real handler registrations the app performs -- the "Apply
# Pipeline" observer in .setup_normalization_handlers() and the "Compare
# Pipelines" observer in .setup_raw_viz_handlers() -- in a minimal server so
# shiny::testServer() can drive them, using the app's own error reporter
# (.show_error_notification) so notifications are exercised for real.
vsn_test_server <- function(notes) {
  function(input, output, session) {
    values <- shiny::reactiveValues(
      raw_data = NULL, normalized_data = NULL,
      pipeline1_data = NULL, pipeline2_data = NULL,
      pipeline1_methods = NULL, pipeline2_methods = NULL,
      current_norm_plot = NULL, current_raw_plot = NULL,
      current_plot = NULL, current_pipeline_plots = NULL,
      plot_history = list()
    )
    add_to_history <- function(...) invisible(NULL)
    LIPIDIFy:::.setup_normalization_handlers(
      input, output, session, values,
      LIPIDIFy:::.show_error_notification, add_to_history
    )
    LIPIDIFy:::.setup_raw_viz_handlers(
      input, output, session, values,
      LIPIDIFy:::.show_error_notification, add_to_history
    )
    notes$values <- values
    invisible(NULL)
  }
}

vsn_capture_notifications <- function(notes) {
  notes$msgs <- character(0)
  testthat::local_mocked_bindings(
    showNotification = function(ui, ...) {
      notes$msgs <- c(notes$msgs, paste(as.character(ui), collapse = " "))
      invisible("mock-id")
    },
    .package = "shiny",
    .env = parent.frame()
  )
}

vsn_raw_data <- function() {
  list(
    numeric_data = vsn_matrix(),
    metadata = data.frame(
      `Sample Name` = paste0("Sample_", seq_len(8)),
      `Sample Group` = rep(c("A", "B"), each = 4),
      check.names = FALSE
    )
  )
}

testthat::test_that("Shiny server dispatch calls VSN when VSN is selected", {
  testthat::skip_if_not_installed("vsn")
  testthat::skip_if_not_installed("shiny")

  called <- new.env(parent = emptyenv())
  called$methods <- character(0)
  testthat::local_mocked_bindings(
    normalize_vsn = function(data, ...) {
      called$methods <- c(called$methods, "VSN")
      data
    },
    normalize_log2median = function(data) {
      called$methods <- c(called$methods, "Log2Median")
      data
    },
    .package = "LIPIDIFy"
  )

  notes <- new.env(parent = emptyenv())
  vsn_capture_notifications(notes)

  shiny::testServer(vsn_test_server(notes), {
    notes$values$raw_data <- vsn_raw_data()
    session$setInputs(
      norm_methods_1 = "VSN", chosen_pipeline = "1",
      apply_normalization = 1
    )
    testthat::expect_identical(called$methods, "VSN")
    testthat::expect_false("Log2Median" %in% called$methods)
    testthat::expect_false(is.null(notes$values$normalized_data))
    testthat::expect_identical(notes$values$normalized_data$methods, "VSN")
    testthat::expect_true(any(grepl(
      "applied successfully: VSN", notes$msgs,
      fixed = TRUE
    )))
  })
})

testthat::test_that("Shiny server dispatch calls Log2Median (not VSN) when Log2Median is selected", {
  testthat::skip_if_not_installed("shiny")

  called <- new.env(parent = emptyenv())
  called$methods <- character(0)
  testthat::local_mocked_bindings(
    normalize_vsn = function(data, ...) {
      called$methods <- c(called$methods, "VSN")
      data
    },
    normalize_log2median = function(data) {
      called$methods <- c(called$methods, "Log2Median")
      data
    },
    .package = "LIPIDIFy"
  )

  notes <- new.env(parent = emptyenv())
  vsn_capture_notifications(notes)

  shiny::testServer(vsn_test_server(notes), {
    notes$values$raw_data <- vsn_raw_data()
    session$setInputs(
      norm_methods_1 = "Log2Median", chosen_pipeline = "1",
      apply_normalization = 1
    )
    testthat::expect_identical(called$methods, "Log2Median")
    testthat::expect_false("VSN" %in% called$methods)
  })
})

testthat::test_that("the app shows the dependency error and never claims VSN succeeded", {
  testthat::skip_if_not_installed("shiny")

  testthat::local_mocked_bindings(
    .vsn_installed = function() FALSE,
    .package = "LIPIDIFy"
  )
  notes <- new.env(parent = emptyenv())
  vsn_capture_notifications(notes)

  shiny::testServer(vsn_test_server(notes), {
    notes$values$raw_data <- vsn_raw_data()
    session$setInputs(
      norm_methods_1 = "VSN", chosen_pipeline = "1",
      apply_normalization = 1
    )

    # No normalized data was committed: nothing silently fell back to
    # Log2Median or any other method.
    testthat::expect_null(notes$values$normalized_data)
    # The dependency error is surfaced to the user verbatim...
    testthat::expect_true(any(grepl(
      "requires the Bioconductor package 'vsn'", notes$msgs,
      fixed = TRUE
    )))
    testthat::expect_true(any(grepl(
      "BiocManager::install('vsn')", notes$msgs,
      fixed = TRUE
    )))
    # ...and success is never reported.
    testthat::expect_false(any(grepl("applied successfully", notes$msgs, fixed = TRUE)))
  })
})

testthat::test_that("pipeline comparison commits nothing when a pipeline fails on missing vsn", {
  testthat::skip_if_not_installed("shiny")

  testthat::local_mocked_bindings(
    .vsn_installed = function() FALSE,
    .package = "LIPIDIFy"
  )
  notes <- new.env(parent = emptyenv())
  vsn_capture_notifications(notes)

  shiny::testServer(vsn_test_server(notes), {
    notes$values$raw_data <- vsn_raw_data()
    session$setInputs(
      norm_methods_1 = "Median", norm_methods_2 = "VSN",
      compare_pipelines = 1
    )
    # Pipeline 1 succeeded on its own, but the run as a whole failed, so
    # neither result is committed and no completion message is shown.
    testthat::expect_null(notes$values$pipeline1_data)
    testthat::expect_null(notes$values$pipeline2_data)
    testthat::expect_false(any(grepl("comparison complete", notes$msgs, fixed = TRUE)))
    testthat::expect_true(any(grepl("BiocManager::install('vsn')", notes$msgs, fixed = TRUE)))
  })
})
