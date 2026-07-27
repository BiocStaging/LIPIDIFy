# =============================================================================
# tests/testthat/test-report-escaping.R
#
# Regression coverage for the report generator's LaTeX-escaping fix.
# Underscores (and other LaTeX-special characters) in dynamically generated
# contrast/lipid/group names -- e.g. "A375_ND_Vehicle" -- previously reached
# raw fig.cap chunk options and Markdown headings unescaped, producing a
# "Missing $ inserted" LaTeX compilation error when the Shiny app's report
# was rendered to PDF (HTML rendering was unaffected).
#
# A second-round failure ("Missing $ inserted. \caption{Expression PC
# 14:1_19:5}") came from a lipid/feature-name caption on an expression plot,
# reached via the *same* extra_plots/embed_plot() code path as the
# originally reported volcano-plot caption -- both are covered explicitly
# below, plus a structural test that fails if any *future* change reintroduces
# a literal (non-variable) fig.cap value anywhere in the generated report.
# =============================================================================

tricky_names <- c(
  "A375_ND_Vehicle",
  "Group_A_vs_Group_B",
  "PC(16:0_18:1)",
  "10%_Treatment",
  "A&B",
  "PC 14:1_19:5"
)

# 1. .escape_report_text (internal) -------------------------------------------
testthat::test_that(".escape_report_text leaves HTML text unchanged", {
  for (nm in tricky_names) {
    testthat::expect_identical(LIPIDIFy:::.escape_report_text(nm, "html"), nm)
  }
  testthat::expect_identical(LIPIDIFy:::.escape_report_text(NULL, "html"), NULL)
})

testthat::test_that(".escape_report_text escapes every required LaTeX-special character", {
  esc <- function(x) LIPIDIFy:::.escape_report_text(x, "pdf")

  testthat::expect_equal(esc("A375_ND_Vehicle"), "A375\\_ND\\_Vehicle")
  testthat::expect_equal(esc("Group_A_vs_Group_B"), "Group\\_A\\_vs\\_Group\\_B")
  testthat::expect_equal(esc("PC(16:0_18:1)"), "PC(16:0\\_18:1)")
  testthat::expect_equal(esc("10%_Treatment"), "10\\%\\_Treatment")
  testthat::expect_equal(esc("A&B"), "A\\&B")
  testthat::expect_equal(esc("PC 14:1_19:5"), "PC 14:1\\_19:5")
  testthat::expect_equal(esc(paste("Expression", "PC 14:1_19:5")), "Expression PC 14:1\\_19:5")

  # every character named in the requirements, in one string
  all_chars <- "back\\slash under_score percent% amp&ersand #hash $dollar {brace} ~tilde ^caret"
  out <- esc(all_chars)
  testthat::expect_true(grepl("\\textbackslash{}", out, fixed = TRUE))
  testthat::expect_true(grepl("\\_", out, fixed = TRUE))
  testthat::expect_true(grepl("\\%", out, fixed = TRUE))
  testthat::expect_true(grepl("\\&", out, fixed = TRUE))
  testthat::expect_true(grepl("\\#", out, fixed = TRUE))
  testthat::expect_true(grepl("\\$", out, fixed = TRUE))
  testthat::expect_true(grepl("\\{", out, fixed = TRUE))
  testthat::expect_true(grepl("\\}", out, fixed = TRUE))
  testthat::expect_true(grepl("\\textasciitilde{}", out, fixed = TRUE))
  testthat::expect_true(grepl("\\textasciicircum{}", out, fixed = TRUE))
})

testthat::test_that(".escape_report_text does not alter the underlying identifier used for lookups", {
  # The fix must only affect display text; contrast names must still be
  # usable, unescaped, as list keys.
  results <- stats::setNames(list(1, 2), tricky_names[1:2])
  for (nm in tricky_names[1:2]) {
    testthat::expect_false(is.null(results[[nm]]))
  }
})

# 2. Build a report Rmd containing every tricky name (including the reported
#    expression-plot lipid caption) via every add_to_history()-fed plot type,
#    and render it -------------------------------------------------------

build_fake_report_inputs <- function(report_dir) {
  make_png <- function(name) {
    fp <- file.path(report_dir, name)
    grDevices::png(fp, width = 200, height = 150)
    graphics::plot(1:5)
    grDevices::dev.off()
    fp
  }

  volcano_png <- make_png("volcano.png")
  expression_png <- make_png("expression.png")

  metadata <- data.frame(
    `Sample Name` = paste0("S", 1:4),
    `Sample Group` = c(
      "A375_ND_Vehicle", "A375_ND_Vehicle",
      "A375_HFD_Vehicle", "A375_HFD_Vehicle"
    ),
    check.names = FALSE
  )
  numeric_data <- matrix(
    stats::rlnorm(4 * 3, 8, 1),
    nrow = 4,
    dimnames = list(
      metadata$`Sample Name`,
      c("PC(16:0_18:1)", "10%_Treatment", "A&B")
    )
  )
  raw_data <- list(
    data = cbind(metadata, numeric_data),
    metadata = metadata,
    numeric_data = numeric_data
  )

  contrast_nm <- "A375_ND_Vehicle - A375_HFD_Vehicle"
  diff_df <- data.frame(
    logFC = c(1.2, -0.8, 0.3),
    AveExpr = c(5, 6, 7),
    P.Value = c(0.001, 0.02, 0.5),
    adj.P.Val = c(0.01, 0.04, 0.6),
    row.names = c("PC(16:0_18:1)", "10%_Treatment", "A&B")
  )
  diff_results <- list(method = "limma", results = stats::setNames(list(diff_df), contrast_nm))

  enrich_df <- data.frame(
    pathway = c("A&B", "PC(16:0_18:1)"),
    NES = c(1.5, -1.1),
    pval = c(0.01, 0.03),
    padj = c(0.05, 0.1)
  )
  enrichment_results <- stats::setNames(
    list(stats::setNames(list(enrich_df), "Group_A_vs_Group_B")),
    contrast_nm
  )

  # extra_plots mirrors exactly what the download_report handler builds from
  # values$plot_history -- one volcano entry (the originally reported
  # caption) and one expression-plot entry using the exact reported lipid
  # name "PC 14:1_19:5" (the second-round failure).
  extra_plots <- list(
    list(
      file = basename(volcano_png),
      label = paste("Volcano", contrast_nm),
      section = "Results"
    ),
    list(
      file = basename(expression_png),
      label = paste("Expression", "PC 14:1_19:5"),
      section = "Lipid Expression"
    )
  )

  list(
    raw_data = raw_data,
    normalized_data = raw_data,
    diff_results = diff_results,
    enrichment_results = enrichment_results,
    extra_plots = extra_plots
  )
}

testthat::test_that("no dynamic fig.cap value is ever inlined as literal chunk-option text", {
  # Structural safeguard: every fig.cap must reference a report_caption_N
  # variable defined in its own setup chunk -- never a literal "..." string,
  # which would be re-parsed as R source and choke on LaTeX escape sequences
  # like "\_" (not a valid R string escape).
  report_dir <- tempfile("lipidify_report_structural_")
  dir.create(report_dir)
  inputs <- build_fake_report_inputs(report_dir)

  for (fmt in c("html", "pdf")) {
    rmd <- LIPIDIFy:::.build_report_rmd_with_plots(
      title              = "Structural Test",
      author             = "Author",
      sections           = c("data_summary", "normalization", "diff_analysis", "enrichment"),
      raw_data           = inputs$raw_data,
      normalized_data    = inputs$normalized_data,
      diff_results       = inputs$diff_results,
      enrichment_results = inputs$enrichment_results,
      output_format      = fmt,
      plot_files         = list(),
      extra_plots        = inputs$extra_plots
    )

    lines <- strsplit(rmd, "\n", fixed = TRUE)[[1]]
    figcap_lines <- grep("fig\\.cap\\s*=", lines, value = TRUE)
    testthat::expect_true(length(figcap_lines) > 0, info = paste("format:", fmt))

    for (l in figcap_lines) {
      # Every fig.cap value must be a bare identifier (report_caption_<n>),
      # never a quoted literal string.
      testthat::expect_false(
        grepl('fig\\.cap\\s*=\\s*"', l),
        info = paste0("Literal fig.cap string found in ", fmt, " output: ", l)
      )
      testthat::expect_true(
        grepl("fig\\.cap\\s*=\\s*report_caption_[0-9]+\\s*[},]", l),
        info = paste0("fig.cap does not reference a report_caption_N variable in ", fmt, " output: ", l)
      )
    }

    # Every referenced report_caption_N variable must have a corresponding
    # assignment chunk earlier in the document.
    var_names <- unique(regmatches(rmd, gregexpr("report_caption_[0-9]+", rmd))[[1]])
    for (v in var_names) {
      testthat::expect_true(
        grepl(paste0(v, ' <- "'), rmd, fixed = TRUE),
        info = paste("Missing assignment chunk for", v)
      )
    }
  }
})

testthat::test_that("generated report escapes all dynamic names (contrast, lipid, group) and renders in both HTML and PDF", {
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")

  report_dir <- tempfile("lipidify_report_")
  dir.create(report_dir)
  inputs <- build_fake_report_inputs(report_dir)

  for (fmt in c("html", "pdf")) {
    if (identical(fmt, "pdf")) {
      testthat::skip_if_not(
        requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex(),
        "tinytex not available"
      )
    }

    rmd <- LIPIDIFy:::.build_report_rmd_with_plots(
      title              = "Test Report: Group_A_vs_Group_B",
      author             = "A&B Author",
      sections           = c("data_summary", "normalization", "diff_analysis", "enrichment"),
      raw_data           = inputs$raw_data,
      normalized_data    = inputs$normalized_data,
      diff_results       = inputs$diff_results,
      enrichment_results = inputs$enrichment_results,
      output_format      = fmt,
      plot_files         = list(),
      extra_plots        = inputs$extra_plots
    )

    # Every tricky name must appear ONLY in its escaped (pdf) or literal
    # (html) form -- never as raw unescaped underscores inside a fig.cap or
    # heading destined for LaTeX.
    if (identical(fmt, "pdf")) {
      testthat::expect_true(grepl("A375\\_ND\\_Vehicle", rmd, fixed = TRUE))
      testthat::expect_true(grepl("PC(16:0\\_18:1)", rmd, fixed = TRUE))
      testthat::expect_true(grepl("10\\%\\_Treatment", rmd, fixed = TRUE))
      testthat::expect_true(grepl("A\\&B", rmd, fixed = TRUE))
      testthat::expect_true(grepl("Expression PC 14:1\\_19:5", rmd, fixed = TRUE))
    } else {
      testthat::expect_true(grepl("A375_ND_Vehicle", rmd, fixed = TRUE))
      testthat::expect_true(grepl("PC(16:0_18:1)", rmd, fixed = TRUE))
      testthat::expect_true(grepl("Expression PC 14:1_19:5", rmd, fixed = TRUE))
    }

    # The YAML title/author must NOT be LaTeX-escaped (that would break YAML
    # parsing outright); pandoc already escapes metadata fields correctly for
    # either output format on its own.
    testthat::expect_true(grepl('title: "Test Report: Group_A_vs_Group_B"', rmd, fixed = TRUE))
    testthat::expect_true(grepl('author: "A&B Author"', rmd, fixed = TRUE))

    rmd_path <- file.path(report_dir, paste0("report_", fmt, ".Rmd"))
    con <- file(rmd_path, "w", encoding = "UTF-8")
    writeLines(rmd, con)
    close(con)

    out_format <- if (identical(fmt, "pdf")) {
      rmarkdown::pdf_document(toc = TRUE, keep_tex = TRUE)
    } else {
      rmarkdown::html_document(toc = TRUE)
    }

    out_file <- testthat::expect_no_error(
      rmarkdown::render(rmd_path,
        output_format = out_format,
        output_dir = report_dir, quiet = TRUE
      )
    )

    testthat::expect_true(file.exists(out_file))
    testthat::expect_gt(file.info(out_file)$size, 0)

    if (identical(fmt, "pdf")) {
      log_file <- file.path(report_dir, paste0("report_", fmt, ".log"))
      if (file.exists(log_file)) {
        log_txt <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
        testthat::expect_false(grepl("Missing $ inserted", log_txt, fixed = TRUE))
        testthat::expect_false(grepl("Undefined control sequence", log_txt, fixed = TRUE))
      }
      tex_file <- file.path(report_dir, paste0("report_", fmt, ".tex"))
      if (file.exists(tex_file)) {
        tex_txt <- paste(readLines(tex_file, warn = FALSE), collapse = "\n")
        testthat::expect_true(grepl("Expression PC 14:1\\_19:5", tex_txt, fixed = TRUE))
      }

      pdftotext_bin <- Sys.which("pdftotext")
      if (nzchar(pdftotext_bin)) {
        txt_out <- tempfile(fileext = ".txt")
        system2(pdftotext_bin, c(shQuote(out_file), shQuote(txt_out)))
        pdf_text <- paste(readLines(txt_out, warn = FALSE), collapse = "\n")
        testthat::expect_true(grepl("PC 14:1_19:5", pdf_text, fixed = TRUE))
        testthat::expect_true(grepl("A375_ND_Vehicle", pdf_text, fixed = TRUE))
      }
    } else {
      html_txt <- paste(readLines(out_file, warn = FALSE), collapse = "\n")
      testthat::expect_true(grepl("PC 14:1_19:5", html_txt, fixed = TRUE))
    }
  }
})
