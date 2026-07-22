# =============================================================================
# tests/testthat/test-shinytest2.R
#
# Browser-driven integration coverage for launch_lipidomics_app(), using a
# real headless Chrome session (via shinytest2/chromote) rather than direct
# function calls. This exercises the app's reactive graph -- including
# server-side observers, dynamically rendered UI (renderUI-based checkbox
# groups and selectize choices, which only materialise once their tab has
# been made visible), and the report-generation download handlers -- in a
# way unit tests calling internal functions directly cannot.
#
# In particular, this test drives the *exact* real-world scenario that
# previously produced "Missing $ inserted" when downloading a PDF report
# with a multi-lipid Expression plot (caption "Expression <lipid name>",
# e.g. "Expression PC 14:1_19:5") through the live app UI: selecting two
# lipids, generating the plot, and downloading both an HTML and a PDF
# report end-to-end.
#
# Skipped entirely (not a failure) when a browser, shinytest2, or tinytex
# is unavailable, or when running on CRAN/Bioconductor build checks where
# a headless browser cannot be assumed.
# =============================================================================

skip_if_no_chrome <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("shinytest2")
  testthat::skip_if_not_installed("chromote")
  chrome_path <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
  testthat::skip_if(is.null(chrome_path) || !nzchar(chrome_path), "Chrome not available")
}

goto_tab <- function(app, tab) {
  app$run_js(sprintf('document.querySelector(\'a[data-value="%s"]\').click();', tab))
  app$wait_for_idle(timeout = 5000)
}

testthat::test_that("launch_lipidomics_app() runs the full analysis pipeline end-to-end in a real browser, including the reported multi-lipid Expression-plot PDF caption bug", {
  skip_if_no_chrome()
  withr_env <- Sys.getenv("NOT_CRAN")
  Sys.setenv(NOT_CRAN = "true")
  on.exit(Sys.setenv(NOT_CRAN = withr_env), add = TRUE)

  app <- shinytest2::AppDriver$new(
    app_dir = launch_lipidomics_app(),
    name = "lipidify-e2e",
    height = 900, width = 1400,
    load_timeout = 30000,
    timeout = 20000
  )
  on.exit(try(app$stop(), silent = TRUE), add = TRUE)

  # ---- Load data and run the full pipeline -------------------------------
  app$click("load_example")
  app$wait_for_idle(timeout = 15000)
  testthat::expect_equal(app$get_value(input = "group_column"), "Sample Group")

  goto_tab(app, "normalization")
  app$click("apply_normalization")
  app$wait_for_idle(timeout = 15000)

  # ---- Multi-lipid expression plot: the exact reported bug's code path ---
  goto_tab(app, "lipid_expression")
  lipid_choices_raw <- app$get_js('
    var el = document.getElementById("selected_lipids");
    Object.keys(el.selectize.options).join("|||");
  ')
  lipid_choices <- strsplit(lipid_choices_raw, "|||", fixed = TRUE)[[1]]
  lipid_choices <- lipid_choices[grepl("^[A-Za-z]+ [0-9]", lipid_choices)]
  testthat::expect_true(length(lipid_choices) >= 2)
  chosen_lipids <- utils::head(lipid_choices, 2)

  app$set_inputs(selected_lipids = chosen_lipids)
  app$wait_for_idle(timeout = 3000)
  app$click("create_expression_plot")
  app$wait_for_idle(timeout = 15000)

  goto_tab(app, "diff_analysis")
  app$click("run_diff_analysis")
  app$wait_for_idle(timeout = 15000)
  testthat::expect_true(nzchar(app$get_value(input = "contrast_select")))

  goto_tab(app, "results_viz")
  app$click("create_viz")
  app$wait_for_idle(timeout = 15000)

  goto_tab(app, "enrichment")
  app$click("run_enrichment")
  app$wait_for_idle(timeout = 20000)

  goto_tab(app, "enrichment_viz")
  app$click("create_enrichment_viz")
  app$wait_for_idle(timeout = 15000)

  # ---- Report generation: HTML and PDF, with a LaTeX-special author name -
  goto_tab(app, "report")
  app$set_inputs(report_author = "A&B_Author")
  app$wait_for_idle(timeout = 2000)

  report_dir <- tempfile("lipidify_e2e_report_")
  dir.create(report_dir)

  html_path <- app$get_download(
    output = "download_report",
    filename = file.path(report_dir, "report.html")
  )
  testthat::expect_true(file.exists(html_path))
  testthat::expect_gt(file.info(html_path)$size, 0)
  html_txt <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl(chosen_lipids[1], html_txt, fixed = TRUE))
  testthat::expect_true(grepl(chosen_lipids[2], html_txt, fixed = TRUE))

  testthat::skip_if_not(
    requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex(),
    "tinytex not available for the PDF portion of this test"
  )

  app$set_inputs(report_format = "pdf")
  app$wait_for_idle(timeout = 2000)
  pdf_path <- app$get_download(
    output = "download_report",
    filename = file.path(report_dir, "report.pdf")
  )
  testthat::expect_true(file.exists(pdf_path))
  testthat::expect_gt(file.info(pdf_path)$size, 10000) # a real rendered PDF, not a tiny error fallback

  pdftotext_bin <- Sys.which("pdftotext")
  if (nzchar(pdftotext_bin)) {
    txt_out <- tempfile(fileext = ".txt")
    system2(pdftotext_bin, c(shQuote(pdf_path), shQuote(txt_out)))
    pdf_text <- paste(readLines(txt_out, warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("Report Error", pdf_text, fixed = TRUE))
    testthat::expect_true(grepl(paste0("Expression ", chosen_lipids[1]), pdf_text, fixed = TRUE))
    testthat::expect_true(grepl(paste0("Expression ", chosen_lipids[2]), pdf_text, fixed = TRUE))
    testthat::expect_true(grepl("A&B_Author", pdf_text, fixed = TRUE))
  }
})
