# =============================================================================
# tests/testthat/test-analysis.R
# Coverage for R/analysis.R: differential analysis, PCA, PLS-DA, enrichment.
# =============================================================================

make_test_data <- function() {
  df <- generate_example_data()
  load_lipidomics_data_from_df(df)
}

# 1. perform_differential_analysis (limma) ----------------------------------
testthat::test_that("perform_differential_analysis (limma) returns expected structure", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  res <- perform_differential_analysis(norm, d$metadata, "Sample Group", method = "limma")

  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("fit", "results", "design", "contrasts", "method") %in% names(res)))
  testthat::expect_equal(res$method, "limma")
  testthat::expect_true(length(res$results) > 0)
  first <- res$results[[1]]
  testthat::expect_true(all(c("logFC", "AveExpr", "P.Value", "adj.P.Val") %in% colnames(first)))
})

testthat::test_that("perform_differential_analysis rejects invalid method", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(
    perform_differential_analysis(norm, d$metadata, "Sample Group", method = "not_a_method"),
    "limma.*edger"
  )
})

testthat::test_that("perform_differential_analysis rejects non-data.frame metadata", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(
    perform_differential_analysis(norm, as.list(d$metadata), "Sample Group"),
    "data.frame"
  )
})

testthat::test_that("perform_differential_analysis rejects missing group_column", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(
    perform_differential_analysis(norm, d$metadata, "Nonexistent Column"),
    "group_column"
  )
})

testthat::test_that("perform_differential_analysis requires contrasts_list with custom design", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  groups <- factor(make.names(d$metadata$`Sample Group`))
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)

  testthat::expect_error(
    perform_differential_analysis(norm, d$metadata, "Sample Group", design = design),
    "contrasts_list"
  )
})

testthat::test_that("perform_differential_analysis accepts a custom design matrix (limma)", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  groups <- factor(make.names(d$metadata$`Sample Group`))
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  rownames(design) <- d$metadata$`Sample Name`

  res <- perform_differential_analysis(
    norm, d$metadata, "Sample Group",
    contrasts_list = c("GroupB - GroupA"),
    method = "limma", design = design
  )
  testthat::expect_equal(names(res$results), "GroupB - GroupA")
  testthat::expect_true(is.null(res$level_mapping))
})

testthat::test_that("perform_differential_analysis (edger) returns expected structure", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  res <- suppressWarnings(
    perform_differential_analysis(norm, d$metadata, "Sample Group", method = "edger")
  )
  testthat::expect_equal(res$method, "edger")
  first <- res$results[[1]]
  testthat::expect_true(all(c("logFC", "AveExpr", "P.Value", "adj.P.Val") %in% colnames(first)))
})

testthat::test_that("perform_differential_analysis (edger) warns about count-data assumption", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_warning(
    perform_differential_analysis(norm, d$metadata, "Sample Group", method = "edger"),
    "integer count data"
  )
})

# 2. perform_pca --------------------------------------------------------------
testthat::test_that("perform_pca returns expected structure", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  res <- perform_pca(norm, d$metadata, "Sample Group")

  testthat::expect_true(all(c("pca_results", "pca_data", "plot", "variance_explained") %in% names(res)))
  testthat::expect_equal(nrow(res$pca_data), nrow(d$metadata))
  testthat::expect_true(all(c("Sample", "PC1", "PC2", "Group") %in% colnames(res$pca_data)))
})

testthat::test_that("perform_pca rejects non-data.frame metadata", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(perform_pca(norm, list(a = 1), "Sample Group"), "data.frame")
})

# 3. perform_plsda ------------------------------------------------------------
testthat::test_that("perform_plsda returns expected structure", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  res <- perform_plsda(norm, d$metadata, "Sample Group")

  testthat::expect_true(all(c("plsda_results", "scores_data", "plot") %in% names(res)))
  testthat::expect_true(all(c("Sample", "Comp1", "Comp2", "Group") %in% colnames(res$scores_data)))
})

testthat::test_that("perform_plsda rejects non-data.frame metadata", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(perform_plsda(norm, list(a = 1), "Sample Group"), "data.frame")
})

# 4. perform_enrichment_analysis ----------------------------------------------
testthat::test_that("perform_enrichment_analysis returns per-contrast, per-category results", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  cls <- classify_lipids(colnames(norm))
  diff <- perform_differential_analysis(norm, d$metadata, "Sample Group", method = "limma")

  enrich <- perform_enrichment_analysis(diff$results, cls, min_set_size = 3, max_set_size = 500)
  testthat::expect_equal(names(enrich), names(diff$results))
  testthat::expect_true("LipidGroup" %in% names(enrich[[1]]))
})

testthat::test_that("perform_enrichment_analysis rejects non-list results_list", {
  cls <- classify_lipids(c("PC 16:0_18:1", "PE 18:0_20:4"))
  testthat::expect_error(perform_enrichment_analysis(data.frame(x = 1), cls), "results_list")
})

testthat::test_that("perform_enrichment_analysis rejects classification_data without Lipid column", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  diff <- perform_differential_analysis(norm, d$metadata, "Sample Group", method = "limma")
  testthat::expect_error(
    perform_enrichment_analysis(diff$results, data.frame(NotLipid = "x")),
    "Lipid"
  )
})

# 5. create_default_contrasts (edge cases not already covered) ---------------
testthat::test_that("create_default_contrasts drops NA levels", {
  ctrs <- create_default_contrasts(c("A", "B", NA))
  testthat::expect_equal(length(ctrs), 1L)
})
