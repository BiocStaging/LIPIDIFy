# =============================================================================
# tests/testthat/test-shiny-server.R
#
# In-process coverage for the Shiny app returned by launch_lipidomics_app().
#
# Unlike test-shinytest2.R, which drives the app in a separate R process
# through a headless browser (and therefore is invisible to covr), these tests
# run the app's real server function with shiny::testServer(). Inputs are set
# with session$setInputs(), which fires the same observers, render functions
# and download handlers the browser would, and the server's `values`
# reactiveValues can be inspected directly. No browser is required, so these
# tests also run during R CMD check and on the Bioconductor build machines.
#
# Note: update*Input() calls do not change inputs inside testServer(), so any
# input that the app would normally populate itself (group columns, contrast
# choices, ...) is set explicitly in the tests below.
# =============================================================================

# Records every notification and modal the server shows, so tests can assert
# on user-facing feedback, including error and warning paths. Also opens a null
# graphics device, so grobs built outside a render function (e.g. for
# downloads) do not leave an Rplots.pdf behind.
capture_ui_events <- function(env = parent.frame()) {
  withr::local_pdf(NULL, .local_envir = env)
  events <- new.env(parent = emptyenv())
  events$notes <- character(0)
  events$types <- character(0)
  events$modals <- character(0)
  testthat::local_mocked_bindings(
    showNotification = function(ui, action = NULL, duration = 5,
                                closeButton = TRUE, id = NULL,
                                type = "default", session = NULL) {
      events$notes <- c(events$notes, paste(as.character(ui), collapse = " "))
      events$types <- c(events$types, type)
      invisible("mock-id")
    },
    removeNotification = function(id, session = NULL) invisible(id),
    showModal = function(ui, session = NULL) {
      events$modals <- c(events$modals, paste(as.character(ui), collapse = " "))
      invisible(NULL)
    },
    .package = "shiny",
    .env = env
  )
  events
}

last_note <- function(events) utils::tail(events$notes, 1)

# Replaces a LIPIDIFy function with a wrapper that calls the real function
# until `failure$message` is set, and then fails with that message. The mock
# must be created at test scope (not inside testServer()), or it outlives the
# test.
mock_failure <- function(name, env = parent.frame()) {
  failure <- new.env(parent = emptyenv())
  failure$message <- NULL
  real_fun <- utils::getFromNamespace(name, "LIPIDIFy")
  mock <- stats::setNames(list(function(...) {
    if (!is.null(failure$message)) stop(failure$message)
    real_fun(...)
  }), name)
  do.call(testthat::local_mocked_bindings, c(mock, .package = "LIPIDIFy", .env = env))
  failure
}

# Loads the example data and applies a Median -> Log2Median pipeline, setting
# the inputs the app's update*Input() calls would otherwise have filled in.
load_and_normalize <- function(session) {
  session$setInputs(load_example = 1)
  session$setInputs(
    group_column = "Sample Group",
    group_col_diff = "Sample Group",
    norm_methods_1 = c("Median", "Log2Median"),
    norm_methods_2 = "TIC",
    chosen_pipeline = "1",
    apply_normalization = 1
  )
}

# Runs limma on two contrasts after load_and_normalize().
run_diff <- function(session) {
  session$setInputs(
    selected_contrasts = c("GroupB - GroupA", "GroupC - GroupA"),
    custom_contrasts_text = "",
    diff_method = "limma",
    run_diff_analysis = 1
  )
  session$setInputs(
    contrast_display = "GroupB - GroupA",
    contrast_select = "GroupB - GroupA"
  )
}

expect_nonempty_file <- function(path) {
  testthat::expect_true(file.exists(path))
  testthat::expect_gt(file.size(path), 0)
}

# 1. App object and UI ---------------------------------------------------------
testthat::test_that("launch_lipidomics_app() returns an app object with the full dashboard UI", {
  app <- launch_lipidomics_app()
  testthat::expect_s3_class(app, "shiny.appobj")

  ui_html <- as.character(LIPIDIFy:::.build_lipidify_ui())
  for (id in c(
    "load_example", "selected_lipids", "run_diff_analysis",
    "run_enrichment", "download_report", "run_batch_correction"
  )) {
    testthat::expect_true(grepl(paste0('id="', id, '"'), ui_html, fixed = TRUE))
  }
})

testthat::test_that("UI helpers produce buttons, checkboxes and download names", {
  tags <- as.character(LIPIDIFy:::.checkbox_group_with_buttons(
    "grp", "Groups", choices = c("A", "B"), selected = "A"
  ))
  testthat::expect_true(grepl("grp_select_all", tags, fixed = TRUE))
  testthat::expect_true(grepl("grp_deselect_all", tags, fixed = TRUE))

  testthat::expect_identical(LIPIDIFy:::.truncate_sheet_name("A / B"), "A___B")
  testthat::expect_identical(nchar(LIPIDIFy:::.truncate_sheet_name(strrep("x", 40))), 31L)
  testthat::expect_match(LIPIDIFy:::.dl_name("raw", "box", ext = "pdf"), "^raw_box_[0-9]{8}\\.pdf$")
})

# 2. Data loading and classification ----------------------------------------
testthat::test_that("example data keeps metadata columns out of the lipid selector", {
  withr::local_seed(1)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    # Reproduces the reviewer's report: with the upload tab's metadata
    # checkboxes at a partial selection, "Tumour ID" and "Weight (mg)" used to
    # end up as lipid features.
    session$setInputs(metadata_cols = c("Sample Name", "Sample Group"))
    session$setInputs(load_example = 1)

    lipids <- colnames(values$raw_data$numeric_data)
    testthat::expect_false(any(c("Tumour ID", "Weight (mg)", "Sample Name") %in% lipids))
    testthat::expect_true(all(c("Tumour ID", "Weight (mg)") %in% names(values$raw_data$metadata)))
    testthat::expect_true(all(grepl("^[A-Za-z]+ [0-9]", lipids)))
    testthat::expect_identical(values$classification$Lipid, lipids)
    testthat::expect_match(last_note(events), "Example data loaded")

    testthat::expect_match(output$data_summary, "Samples: 20")
    testthat::expect_match(output$data_summary, "GroupA, GroupB, GroupC, GroupD")
    testthat::expect_true(nzchar(output$data_preview))
    testthat::expect_true(nzchar(output$classification_table))
    cls_file <- output$download_classification
    testthat::expect_identical(utils::read.csv(cls_file, check.names = FALSE)$Lipid, lipids)
  })
})

testthat::test_that("the upload tab defaults to all example metadata columns", {
  ui_html <- as.character(LIPIDIFy:::.ui_tab_upload())
  for (col in c("Tumour ID", "Weight (mg)")) {
    testthat::expect_match(
      ui_html,
      paste0('value="', col, '" checked="checked"'),
      fixed = TRUE
    )
  }
})

testthat::test_that("uploaded CSV files load, and unreadable files report an error", {
  withr::local_seed(2)
  csv <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(generate_example_data(), csv, row.names = FALSE)
  events <- capture_ui_events()

  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(
      metadata_cols = LIPIDIFy:::.example_metadata_cols,
      file = list(datapath = csv, name = "example.csv")
    )
    session$setInputs(load_data = 1)
    testthat::expect_identical(nrow(values$raw_data$numeric_data), 20L)
    testthat::expect_false("Weight (mg)" %in% colnames(values$raw_data$numeric_data))
    testthat::expect_match(last_note(events), "Data loaded successfully")

    session$setInputs(file = list(datapath = tempfile(fileext = ".csv"), name = "missing.csv"))
    session$setInputs(load_data = 2)
    testthat::expect_match(last_note(events), "Error loading data")
  })
})

testthat::test_that("custom classification can be loaded, rejected and reset", {
  withr::local_seed(3)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    lipids <- colnames(values$raw_data$numeric_data)

    good <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(
      Lipid = lipids,
      MyClass = sub(" .*", "", lipids)
    ), good, row.names = FALSE)
    session$setInputs(custom_classification_file = list(datapath = good))
    session$setInputs(load_custom_classification = 1)
    testthat::expect_true(values$use_custom_classification)
    testthat::expect_true("MyClass" %in% colnames(values$classification_data))
    testthat::expect_true(nzchar(output$classification_table))
    cls <- utils::read.csv(output$download_classification)
    testthat::expect_true("MyClass" %in% colnames(cls))

    bad <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(Name = "x", Class = "y"), bad, row.names = FALSE)
    session$setInputs(custom_classification_file = list(datapath = bad))
    session$setInputs(load_custom_classification = 2)
    testthat::expect_match(last_note(events), "Loading custom classification failed")

    session$setInputs(reset_classification = 1)
    testthat::expect_false(values$use_custom_classification)
    testthat::expect_null(values$custom_classification)
    testthat::expect_match(last_note(events), "Automatic classification restored")
  })
})

testthat::test_that("select-all / deselect-all buttons and help dialogs respond", {
  withr::local_seed(4)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1, group_column = "Sample Group")
    ids <- c(
      "metadata_cols", "norm_methods_1", "norm_methods_2", "groups_included",
      "selected_samples", "selected_expression_groups", "report_sections",
      "raw_filter_groups", "norm_filter_groups"
    )
    for (id in ids) {
      do.call(session$setInputs, stats::setNames(list(1), paste0(id, "_select_all")))
      do.call(session$setInputs, stats::setNames(list(1), paste0(id, "_deselect_all")))
    }
    testthat::expect_length(events$notes[events$types == "error"], 0L)

    info_ids <- c(
      "info_classification", "info_norm_methods", "info_ellipses",
      "info_methods", "info_contrasts", "info_custom_sets",
      "info_imputation", "info_batch"
    )
    for (id in info_ids) {
      do.call(session$setInputs, stats::setNames(list(1), id))
    }
    testthat::expect_length(events$modals, length(info_ids))
    testthat::expect_true(any(grepl("Quantile Normalization", events$modals)))
  })
})

# 3. Raw data visualization and pipeline comparison ---------------------------
testthat::test_that("raw data plots render for every plot type and download", {
  withr::local_seed(5)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    session$setInputs(
      raw_filter_groups = c("GroupA", "GroupB", "GroupC"),
      raw_color_by_group = TRUE, view_mode = "sample",
      limit_top_samples = TRUE, top_n_samples = 10,
      top_n_lipids = 20, raw_heatmap_top_n = 20,
      raw_group_column = "Sample Group",
      raw_groups_included = c("GroupA", "GroupB"),
      raw_ellipse_type = "confidence", raw_show_sample_labels = TRUE,
      img_format_raw = "png"
    )

    for (pt in c("boxplot", "violin", "density")) {
      for (vm in c("sample", "lipid")) {
        n_hist <- length(values$plot_history)
        session$setInputs(plot_type = pt, view_mode = vm)
        session$setInputs(create_raw_plot = sample.int(1e6, 1))
        testthat::expect_length(values$plot_history, n_hist + 1L)
        testthat::expect_true(nzchar(output$raw_plot))
      }
    }
    expect_nonempty_file(output$download_raw_plot)

    session$setInputs(plot_type = "heatmap", create_raw_plot = 100)
    testthat::expect_s3_class(values$current_raw_plot, "pheatmap")
    testthat::expect_true(inherits(output$raw_heatmap_plot, "list"))
    session$setInputs(img_format_raw = "pdf")
    expect_nonempty_file(output$download_raw_plot)

    session$setInputs(plot_type = "pca", create_raw_plot = 101)
    testthat::expect_true(any(grepl("Raw PCA", vapply(values$plot_history, `[[`, "", "label"))))
    session$setInputs(plot_type = "plsda", create_raw_plot = 102)
    testthat::expect_true(any(grepl("Raw PLS-DA", vapply(values$plot_history, `[[`, "", "label"))))
    testthat::expect_true(nzchar(output$raw_plot))
    testthat::expect_length(events$notes[events$types == "error"], 0L)

    session$setInputs(raw_group_column = "No Such Column", create_raw_plot = 103)
    testthat::expect_match(last_note(events), "was not found in the metadata")
  })
})

testthat::test_that("raw PCA warns about, and raw PLS-DA refuses, missing values", {
  withr::local_seed(6)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    rd <- values$raw_data
    rd$numeric_data[1:3, 1] <- NA
    values$raw_data <- rd
    session$setInputs(
      raw_group_column = "Sample Group", raw_groups_included = NULL,
      raw_ellipse_type = "none", plot_type = "plsda", create_raw_plot = 1
    )
    testthat::expect_match(last_note(events), "PLS-DA cannot be computed on raw data")
    testthat::expect_null(values$current_raw_plot)

    testthat::expect_warning(
      session$setInputs(plot_type = "pca", create_raw_plot = 2),
      "imputed by the mean"
    )
    testthat::expect_true(any(grepl("missing value\\(s\\) were replaced", events$notes)))
    testthat::expect_false(is.null(values$current_raw_plot))
  })
})

testthat::test_that("normalization pipelines are compared, plotted and downloaded", {
  withr::local_seed(7)
  events <- capture_ui_events()
  norm_failure <- mock_failure("apply_normalizations")
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    session$setInputs(
      norm_methods_1 = c("Median", "Log2Median"), norm_methods_2 = "Quantile",
      norm_compare_plot_type = "boxplot", norm_color_by_group = TRUE,
      norm_compare_view_mode = "sample", norm_compare_limit_top = TRUE,
      norm_compare_top_n = 10, img_format_pipeline = "png"
    )
    session$setInputs(compare_pipelines = 1)
    testthat::expect_match(last_note(events), "Pipeline comparison complete")
    testthat::expect_identical(values$pipeline2_methods, "Quantile")

    output$pipeline_comparison
    testthat::expect_named(values$current_pipeline_plots, c("p1", "p2"))
    expect_nonempty_file(output$download_pipeline_comparison)
    session$setInputs(img_format_pipeline = "pdf", norm_compare_limit_top = FALSE)
    expect_nonempty_file(output$download_pipeline_comparison)

    # if either pipeline fails, neither is committed
    norm_failure$message <- "simulated failure"
    session$setInputs(norm_methods_2 = "TIC", compare_pipelines = 2)
    testthat::expect_match(last_note(events), "Pipeline comparison failed: simulated failure")
    testthat::expect_identical(values$pipeline2_methods, "Quantile")
  })
})

# 4. Normalization, preprocessing and normalized plots ------------------------
testthat::test_that("normalized data plots render for every plot type and download", {
  withr::local_seed(8)
  events <- capture_ui_events()
  norm_failure <- mock_failure("apply_normalizations")
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    testthat::expect_identical(values$normalized_data$methods, c("Median", "Log2Median"))
    testthat::expect_match(last_note(events), "Median -> Log2Median")

    session$setInputs(
      norm_filter_groups = c("GroupA", "GroupB"),
      norm_color_by_group_viz = TRUE, norm_limit_top_samples = TRUE,
      norm_top_n_samples = 5, norm_heatmap_top_n = 20,
      ellipse_type = "visual", show_sample_labels = TRUE,
      img_format_norm = "png"
    )
    for (pt in c("boxplot", "violin", "density", "pca", "plsda")) {
      session$setInputs(norm_plot_type = pt)
      session$setInputs(create_norm_plot = sample.int(1e6, 1))
      testthat::expect_true(nzchar(output$norm_plot))
    }
    expect_nonempty_file(output$download_norm_plot)

    # groups_included is used when no norm_filter_groups are selected
    session$setInputs(
      norm_filter_groups = character(0),
      groups_included = c("GroupC", "GroupD"),
      norm_plot_type = "heatmap", create_norm_plot = 100
    )
    testthat::expect_s3_class(values$current_norm_plot, "pheatmap")
    output$norm_heatmap_plot
    session$setInputs(img_format_norm = "pdf")
    expect_nonempty_file(output$download_norm_plot)
    testthat::expect_length(events$notes[events$types == "error"], 0L)

    # the second pipeline can be chosen, and failures are reported
    session$setInputs(chosen_pipeline = "2", apply_normalization = 2)
    testthat::expect_identical(values$normalized_data$methods, "TIC")
    norm_failure$message <- "simulated failure"
    session$setInputs(chosen_pipeline = "1", apply_normalization = 3)
    testthat::expect_match(last_note(events), "Applying normalization failed: simulated failure")
    testthat::expect_identical(values$normalized_data$methods, "TIC")
  })
})

testthat::test_that("imputation is applied, summarised and reset", {
  withr::local_seed(9)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    testthat::expect_match(output$missing_value_summary, "Apply normalisation first")
    session$setInputs(imputation_method = "half_min", imputation_k = 5, run_imputation = 1)
    testthat::expect_match(last_note(events), "apply normalisation first")
    session$setInputs(reset_imputation = 1)
    testthat::expect_match(last_note(events), "No imputation to reset")

    load_and_normalize(session)
    nd <- values$normalized_data
    nd$numeric_data[1:2, 1:2] <- NA
    values$normalized_data <- nd
    session$flushReact()
    testthat::expect_match(output$missing_value_summary, "Missing \\(NA\\)  : 4")
    testthat::expect_match(output$preprocessing_status, "Imputation       : Not applied")

    session$setInputs(imputation_method = "knn", imputation_k = 3, run_imputation = 2)
    testthat::expect_false(anyNA(values$normalized_data$numeric_data))
    testthat::expect_match(output$missing_value_summary, "Imputation    : Applied")
    testthat::expect_match(output$preprocessing_status, "Imputation       : Applied")

    session$setInputs(reset_imputation = 2)
    testthat::expect_true(anyNA(values$normalized_data$numeric_data))
    testthat::expect_null(values$pre_impute_data)

    session$setInputs(imputation_method = "not_a_method", run_imputation = 3)
    testthat::expect_match(last_note(events), "Imputation failed")
  })
})

testthat::test_that("batch correction validates its inputs, runs and resets", {
  withr::local_seed(10)
  events <- capture_ui_events()
  batch_failure <- mock_failure("correct_batch_effects")
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(batch_method = "limma", run_batch_correction = 1)
    testthat::expect_match(last_note(events), "apply normalisation first")

    load_and_normalize(session)
    session$setInputs(batch_column = "Not a column", run_batch_correction = 2)
    testthat::expect_match(last_note(events), "No valid batch column selected")

    md <- values$normalized_data$metadata
    md$Constant <- "run1"
    md$Batch <- rep(c("b1", "b2"), times = 10)
    nd <- values$normalized_data
    nd$metadata <- md
    values$normalized_data <- nd

    session$setInputs(batch_column = "Constant", run_batch_correction = 3)
    testthat::expect_match(last_note(events), "has only one unique value")
    session$setInputs(batch_column = "Sample Name", run_batch_correction = 4)
    testthat::expect_match(last_note(events), "Every sample has a unique value")

    before <- values$normalized_data$numeric_data
    # Batch alternates within every group, so the design is not confounded
    # and the groups stay protected (no warning).
    testthat::expect_no_warning(session$setInputs(
      batch_column = "Batch", batch_group_column = "Sample Group",
      run_batch_correction = 5
    ))
    testthat::expect_match(last_note(events), "Batch correction applied using limma")
    testthat::expect_false(isTRUE(all.equal(before, values$normalized_data$numeric_data)))
    testthat::expect_match(output$preprocessing_status, "Batch correction : Applied")

    session$setInputs(reset_batch_correction = 1)
    testthat::expect_equal(values$normalized_data$numeric_data, before)
    testthat::expect_match(last_note(events), "Batch correction reset")

    # batch column == group column: warned, but still attempted
    testthat::expect_warning(
      session$setInputs(
        batch_column = "Sample Group", batch_group_column = "Sample Group",
        run_batch_correction = 6
      ),
      "confounded"
    )
    testthat::expect_true(any(grepl("batch column and group column are the same", events$notes)))
    session$setInputs(reset_batch_correction = 2)

    # a failing correction is rolled back
    batch_failure$message <- "simulated failure"
    session$setInputs(batch_column = "Batch", batch_group_column = "Sample Group", run_batch_correction = 7)
    testthat::expect_match(last_note(events), "Batch correction failed: simulated failure")
    testthat::expect_equal(values$normalized_data$numeric_data, before)
  })
})

testthat::test_that("batch correction reports a ComBat -> limma fallback", {
  withr::local_seed(11)
  events <- capture_ui_events()
  testthat::local_mocked_bindings(
    correct_batch_effects = function(data_matrix, ...) {
      structure(data_matrix, method_used = "limma")
    },
    .package = "LIPIDIFy"
  )
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    nd <- values$normalized_data
    nd$metadata$Batch <- rep(c("b1", "b2"), times = 10)
    values$normalized_data <- nd
    session$setInputs(
      batch_column = "Batch", batch_group_column = "Sample Group",
      batch_method = "combat", run_batch_correction = 1
    )
    testthat::expect_match(last_note(events), "ComBat requires the 'sva' package")
  })
})

# 5. Lipid expression plots ---------------------------------------------------
testthat::test_that("lipid expression plots render by sample and by group, and download", {
  withr::local_seed(12)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    lipids <- utils::head(colnames(values$raw_data$numeric_data), 2)

    session$setInputs(
      selected_lipids = lipids,
      expression_selection_mode = "samples",
      selected_samples = paste0("Sample_", 1:10),
      expression_data_type = "normalized",
      img_format_expression = "png"
    )
    session$setInputs(create_expression_plot = 1)
    testthat::expect_named(values$current_expression_plots, lipids)
    testthat::expect_true(grepl(make.names(lipids[1]), output$expression_plots_ui$html))
    output[[paste0("expr_plot_", make.names(lipids[1]))]]
    testthat::expect_match(output$plot_history_status, "Expression")
    expect_nonempty_file(output$download_expression_plots)

    session$setInputs(
      expression_selection_mode = "groups",
      selected_expression_groups = c("GroupA", "GroupB"),
      expression_data_type = "raw",
      img_format_expression = "pdf",
      create_expression_plot = 2
    )
    expect_nonempty_file(output$download_expression_plots)

    # a single ggplot (rather than a list) is shown and saved on its own
    values$current_expression_plots <- values$current_expression_plots[[1]]
    session$flushReact()
    testthat::expect_true(grepl("expr_plot_single", output$expression_plots_ui$html))
    output$expr_plot_single
    expect_nonempty_file(output$download_expression_plots)
    testthat::expect_length(events$notes[events$types == "error"], 0L)

    session$setInputs(selected_lipids = "Not a lipid", create_expression_plot = 3)
    testthat::expect_match(last_note(events), "Creating expression plot failed")
  })
})

# 6. Differential analysis and results ----------------------------------------
testthat::test_that("differential analysis runs, summarises and exports results", {
  withr::local_seed(13)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    testthat::expect_match(output$diff_summary, "No results yet")
    load_and_normalize(session)

    testthat::expect_length(values$available_contrasts, 6L)
    testthat::expect_true(grepl("selected_contrasts", output$contrast_selection_ui$html))
    session$setInputs(selected_contrasts_select_all = 1, selected_contrasts_deselect_all = 1)

    session$setInputs(selected_contrasts = NULL, custom_contrasts_text = "  ", run_diff_analysis = 1)
    testthat::expect_match(last_note(events), "Select at least one contrast")
    testthat::expect_null(values$diff_results)

    session$setInputs(
      selected_contrasts = "GroupB - GroupA",
      custom_contrasts_text = "# a comment\n(GroupC + GroupD)/2 - GroupA\n",
      diff_method = "limma", run_diff_analysis = 2
    )
    testthat::expect_setequal(
      names(values$diff_results$results),
      c("GroupB - GroupA", "(GroupC + GroupD)/2 - GroupA")
    )
    testthat::expect_match(output$diff_summary, "=== GroupB - GroupA ===")
    testthat::expect_match(output$report_status, "Diff. analysis done:     TRUE")

    session$setInputs(contrast_display = "GroupB - GroupA")
    testthat::expect_true(nzchar(output$diff_results_table))
    session$setInputs(contrast_display = "missing contrast")
    testthat::expect_true(nzchar(output$diff_results_table))

    session$setInputs(contrast_display = "GroupB - GroupA")
    tbl <- utils::read.csv(output$download_current_table, check.names = FALSE)
    testthat::expect_identical(colnames(tbl)[1], "Lipid")
    testthat::expect_identical(nrow(tbl), ncol(values$normalized_data$numeric_data))

    xlsx <- output$download_results
    testthat::expect_setequal(
      openxlsx::getSheetNames(xlsx),
      vapply(names(values$diff_results$results), LIPIDIFy:::.truncate_sheet_name, "")
    )

    session$setInputs(selected_contrasts = "Nope - Nada", custom_contrasts_text = "", run_diff_analysis = 3)
    testthat::expect_match(last_note(events), "Differential analysis failed")
  })
})

testthat::test_that("edgeR differential analysis runs from the app", {
  withr::local_seed(14)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    # edgeR needs non-negative, count-like data, so use a TIC pipeline
    session$setInputs(
      group_col_diff = "Sample Group", norm_methods_1 = "TIC",
      chosen_pipeline = "1", apply_normalization = 1
    )
    testthat::expect_warning(
      session$setInputs(
        selected_contrasts = "GroupB - GroupA", custom_contrasts_text = "",
        diff_method = "edger", run_diff_analysis = 1
      ),
      "edgeR is designed for integer count data"
    )
    testthat::expect_match(last_note(events), "completed using edger")
    testthat::expect_identical(names(values$diff_results$results), "GroupB - GroupA")
  })
})

testthat::test_that("volcano plots and heatmaps are created from results and downloaded", {
  withr::local_seed(15)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    run_diff(session)

    session$setInputs(
      viz_type = "volcano", logfc_threshold = 0.5, pval_threshold = 0.05,
      top_labels = 5, color_by_class = TRUE, color_column = "LipidGroup",
      img_format_results = "png", create_viz = 1
    )
    testthat::expect_true(inherits(values$current_results_plot, "ggplot"))
    output$results_plot
    expect_nonempty_file(output$download_viz)
    session$setInputs(img_format_results = "pdf")
    expect_nonempty_file(output$download_viz)

    # heatmap of significant lipids, restricted to the two contrast groups
    res <- values$diff_results$results[["GroupB - GroupA"]]
    res$adj.P.Val[1:5] <- 0.001
    res$logFC[1:5] <- 3
    dr <- values$diff_results
    dr$results[["GroupB - GroupA"]] <- res
    values$diff_results <- dr
    session$setInputs(viz_type = "heatmap", heatmap_top_n = 20, create_viz = 2)
    testthat::expect_s3_class(values$current_results_plot, "pheatmap")
    output$results_plot
    expect_nonempty_file(output$download_viz)
    session$setInputs(img_format_results = "png")
    expect_nonempty_file(output$download_viz)

    # significant names that are absent from the data matrix
    rownames(res)[1:5] <- paste0("Unknown_", 1:5)
    dr$results[["GroupB - GroupA"]] <- res
    values$diff_results <- dr
    session$setInputs(create_viz = 3)
    output$results_plot

    # no significant features
    res$adj.P.Val <- 1
    dr$results[["GroupB - GroupA"]] <- res
    values$diff_results <- dr
    session$setInputs(create_viz = 4)
    output$results_plot
    testthat::expect_length(events$notes[events$types == "error"], 0L)
  })
})

# 7. Enrichment ---------------------------------------------------------------
testthat::test_that("enrichment analysis runs, visualises and exports results", {
  withr::local_seed(16)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    run_diff(session)

    lipids <- colnames(values$raw_data$numeric_data)
    sets_csv <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(
      Lipid = lipids,
      Set_Name = ifelse(seq_along(lipids) %% 2 == 0, "Even set", "Odd set")
    ), sets_csv, row.names = FALSE)
    session$setInputs(custom_enrichment_file = list(datapath = sets_csv))
    session$setInputs(load_custom_sets = 1)
    testthat::expect_match(last_note(events), "Custom enrichment sets loaded")

    session$setInputs(min_set_size = 3, max_set_size = 500, run_enrichment = 1)
    testthat::expect_match(last_note(events), "Enrichment analysis completed")
    contrast <- names(values$enrichment_results)[1]
    types <- names(values$enrichment_results[[contrast]])
    testthat::expect_true("LipidGroup" %in% types)

    session$setInputs(enrichment_contrast = contrast, enrichment_type = "LipidGroup")
    testthat::expect_true(nzchar(output$enrichment_results_table))
    enr <- utils::read.csv(output$download_enrichment)
    testthat::expect_true("pathway" %in% colnames(enr))
    session$setInputs(enrichment_type = "NotAType")
    testthat::expect_true(nzchar(output$enrichment_results_table))

    session$setInputs(
      enrichment_viz_contrast = contrast, enrichment_viz_type = "LipidGroup",
      max_pathways = 10, img_format_enrich = "png"
    )
    for (pt in c("dotplot", "barplot")) {
      session$setInputs(enrichment_plot_type = pt, create_enrichment_viz = sample.int(1e6, 1))
      testthat::expect_true(inherits(values$current_enrichment_plot, "ggplot"))
      output$enrichment_plot
    }
    expect_nonempty_file(output$download_enrichment_viz)

    # an empty result table produces a placeholder plot
    er <- values$enrichment_results
    er[[contrast]]$LipidGroup <- er[[contrast]]$LipidGroup[0, ]
    values$enrichment_results <- er
    session$setInputs(create_enrichment_viz = 1e7)
    output$enrichment_plot
    testthat::expect_true(nzchar(output$enrichment_results_table))
    session$setInputs(img_format_enrich = "pdf")
    expect_nonempty_file(output$download_enrichment_viz)
    testthat::expect_length(events$notes[events$types == "error"], 0L)

    bad <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(Lipid = "x"), bad, row.names = FALSE)
    session$setInputs(custom_enrichment_file = list(datapath = bad))
    session$setInputs(load_custom_sets = 2)
    testthat::expect_match(last_note(events), "Loading custom enrichment sets failed")
  })
})

testthat::test_that("enrichment failures are reported, with a hint for fgsea errors", {
  withr::local_seed(17)
  events <- capture_ui_events()
  enrichment_failure <- mock_failure("perform_enrichment_analysis")
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    run_diff(session)
    session$setInputs(min_set_size = 3, max_set_size = 500)

    enrichment_failure$message <- "fgsea is corrupt"
    session$setInputs(run_enrichment = 1)
    testthat::expect_match(last_note(events), "try reinstalling")

    enrichment_failure$message <- "something else"
    session$setInputs(run_enrichment = 2)
    testthat::expect_match(last_note(events), "Enrichment analysis failed: something else")
  })
})

# 8. Report -------------------------------------------------------------------
testthat::test_that("the HTML report includes every analysis step and sessionInfo()", {
  testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  withr::local_seed(18)
  events <- capture_ui_events()
  shiny::testServer(launch_lipidomics_app(), {
    load_and_normalize(session)
    session$setInputs(
      norm_compare_plot_type = "boxplot", norm_color_by_group = FALSE,
      compare_pipelines = 1
    )
    output$pipeline_comparison
    session$setInputs(
      plot_type = "heatmap", raw_filter_groups = NULL,
      raw_heatmap_top_n = 20, create_raw_plot = 1
    )
    session$setInputs(
      norm_plot_type = "boxplot", norm_filter_groups = NULL,
      norm_limit_top_samples = FALSE, create_norm_plot = 1
    )
    run_diff(session)
    session$setInputs(
      viz_type = "volcano", logfc_threshold = 1, pval_threshold = 0.05,
      top_labels = 5, color_by_class = FALSE, create_viz = 1
    )
    session$setInputs(min_set_size = 3, max_set_size = 500, run_enrichment = 1)
    contrast <- names(values$enrichment_results)[1]
    session$setInputs(
      enrichment_viz_contrast = contrast, enrichment_viz_type = "LipidGroup",
      enrichment_plot_type = "barplot", max_pathways = 10,
      create_enrichment_viz = 1
    )

    session$setInputs(
      report_title = "Test \"report\"", report_author = "Tester",
      report_sections = c("data_summary", "normalization", "diff_analysis", "enrichment"),
      report_format = "html"
    )
    report <- output$download_report
    testthat::expect_match(basename(report), "\\.html$")
    html <- paste(readLines(report, warn = FALSE), collapse = "\n")
    testthat::expect_match(html, "Session Information")
    testthat::expect_match(html, "R version", fixed = TRUE)
    testthat::expect_match(html, "attached base packages")
    testthat::expect_match(html, "Additional Visualisations")
    testthat::expect_match(last_note(events), "Report generated successfully")
  })
})

testthat::test_that("report generation errors are written to the download instead of crashing", {
  events <- capture_ui_events()
  testthat::local_mocked_bindings(
    .build_report_rmd_with_plots = function(...) stop("simulated report failure"),
    .package = "LIPIDIFy"
  )
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(
      report_title = "T", report_author = "A",
      report_sections = "data_summary", report_format = "pdf"
    )
    report <- output$download_report
    testthat::expect_match(basename(report), "\\.pdf$")
    testthat::expect_match(paste(readLines(report), collapse = "\n"), "simulated report failure")
    testthat::expect_match(last_note(events), "Report failed")
  })
})

testthat::test_that("the report falls back to HTML when PDF rendering fails", {
  testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  events <- capture_ui_events()
  # Simulates a machine without LaTeX: PDF rendering fails, and the HTML
  # fallback render is recorded rather than run (real HTML rendering is
  # covered above).
  rendered <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    render = function(input, output_format = NULL, output_dir = NULL, ...) {
      if (identical(output_format$pandoc$to, "latex")) stop("no LaTeX here")
      rendered$rmd <- paste(readLines(input), collapse = "\n")
      out <- file.path(output_dir, "report.html")
      writeLines("<html>fallback</html>", out)
      out
    },
    .package = "rmarkdown"
  )
  withr::local_seed(19)
  shiny::testServer(launch_lipidomics_app(), {
    session$setInputs(load_example = 1)
    session$setInputs(
      report_title = "T", report_author = "A",
      report_sections = c("data_summary", "normalization", "diff_analysis", "enrichment"),
      report_format = "pdf"
    )
    report <- output$download_report
    testthat::expect_true(any(grepl("PDF failed - generating HTML instead", events$notes)))
    testthat::expect_identical(readLines(report), "<html>fallback</html>")
    testthat::expect_match(rendered$rmd, "Normalisation has not been applied")
    testthat::expect_match(rendered$rmd, "Differential analysis has not been run")
    testthat::expect_match(rendered$rmd, "Enrichment analysis has not been run")
  })
})

testthat::test_that("the report Rmd contains a sessionInfo() chunk for both formats", {
  for (fmt in c("html", "pdf")) {
    rmd <- LIPIDIFy:::.build_report_rmd_with_plots(
      title = "T", author = "A", sections = character(0),
      raw_data = NULL, normalized_data = NULL,
      diff_results = NULL, enrichment_results = NULL,
      output_format = fmt
    )
    testthat::expect_match(rmd, "```{r session-info", fixed = TRUE)
    testthat::expect_match(rmd, "\nsessionInfo()\n", fixed = TRUE)
  }
})

# 9. Internal plot saver ------------------------------------------------------
testthat::test_that(".save_plot writes ggplot, pheatmap and grob objects, and rejects others", {
  withr::local_seed(20)
  withr::local_pdf(NULL)
  mat <- matrix(stats::rnorm(60), nrow = 10, dimnames = list(paste0("L", 1:10), paste0("S", 1:6)))
  hm <- pheatmap::pheatmap(mat, silent = TRUE)
  gg <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  grob <- gridExtra::arrangeGrob(gg, gg, ncol = 1)

  for (fmt in c("png", "pdf")) {
    for (obj in list(hm, gg, grob)) {
      f <- withr::local_tempfile(fileext = paste0(".", fmt))
      LIPIDIFy:::.save_plot(obj, f, fmt)
      expect_nonempty_file(f)
    }
  }
  testthat::expect_error(
    LIPIDIFy:::.save_plot(list(), tempfile(), "png"),
    "Cannot save object of class"
  )
})
