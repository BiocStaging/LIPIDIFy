# =============================================================================
# tests/testthat/test-visualization.R
# Coverage for R/visualization.R.
# =============================================================================

make_test_data <- function() {
  df <- generate_example_data()
  load_lipidomics_data_from_df(df)
}

# 1. visualize_raw_data -------------------------------------------------------
testthat::test_that("visualize_raw_data returns a ggplot for each plot_type", {
  d <- make_test_data()
  dl <- list(numeric_data = d$numeric_data)
  for (pt in c("boxplot", "density", "histogram")) {
    p <- visualize_raw_data(dl, pt)
    testthat::expect_s3_class(p, "ggplot")
  }
})

testthat::test_that("visualize_raw_data errors on unknown plot_type", {
  d <- make_test_data()
  dl <- list(numeric_data = d$numeric_data)
  testthat::expect_error(visualize_raw_data(dl, "not_a_type"))
})

testthat::test_that("visualize_raw_data rejects a data_list without numeric_data", {
  testthat::expect_error(visualize_raw_data(list(metadata = data.frame())), "numeric_data")
  testthat::expect_error(visualize_raw_data("not a list"), "numeric_data")
})

# 2. visualize_raw_data_improved ----------------------------------------------
testthat::test_that("visualize_raw_data_improved handles sample and lipid view modes", {
  d <- make_test_data()
  dl <- list(numeric_data = d$numeric_data)

  p_sample <- visualize_raw_data_improved(dl, "boxplot", "sample", metadata = d$metadata)
  testthat::expect_s3_class(p_sample, "ggplot")

  p_lipid <- visualize_raw_data_improved(dl, "boxplot", "lipid", top_n = 10)
  testthat::expect_s3_class(p_lipid, "ggplot")
})

testthat::test_that("visualize_raw_data_improved returns placeholder plot for missing numeric_data", {
  p <- visualize_raw_data_improved(list(numeric_data = NULL))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("visualize_raw_data_improved handles unknown view_mode gracefully", {
  d <- make_test_data()
  dl <- list(numeric_data = d$numeric_data)
  p <- visualize_raw_data_improved(dl, "boxplot", "not_a_mode")
  testthat::expect_s3_class(p, "ggplot")
})

# 3. create_pipeline_plot -------------------------------------------------
testthat::test_that("create_pipeline_plot works for boxplot/violin/density and both view modes", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))

  for (pt in c("boxplot", "violin", "density")) {
    p <- create_pipeline_plot(norm, plot_type = pt, metadata = d$metadata)
    testthat::expect_s3_class(p, "ggplot")
  }

  p_lipid <- create_pipeline_plot(norm, plot_type = "boxplot", view_mode = "lipid", top_n = 10)
  testthat::expect_s3_class(p_lipid, "ggplot")
})

testthat::test_that("create_pipeline_plot rejects non-numeric data_matrix", {
  testthat::expect_error(create_pipeline_plot(data.frame(a = "x", b = "y")), "numeric")
})

testthat::test_that("create_pipeline_plot works with metadata lacking the group column", {
  d <- make_test_data()
  p <- create_pipeline_plot(d$numeric_data, metadata = data.frame(Other = seq_len(nrow(d$numeric_data))))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("create_pipeline_plot rejects non-data.frame metadata", {
  d <- make_test_data()
  testthat::expect_error(create_pipeline_plot(d$numeric_data, metadata = list(a = 1)), "data.frame")
})

# 4. create_heatmap_robust -----------------------------------------------
testthat::test_that("create_heatmap_robust returns a pheatmap object", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  hm <- create_heatmap_robust(t(norm), d$metadata, "Sample Group", top_n = 10)
  testthat::expect_s3_class(hm, "pheatmap")
})

testthat::test_that("create_heatmap_robust returns an error plot for bad input", {
  p <- create_heatmap_robust(matrix(nrow = 0, ncol = 0), data.frame(x = 1))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("create_heatmap_robust returns an error plot for invalid metadata", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  p <- create_heatmap_robust(t(norm), list(not = "a data.frame"), "Sample Group")
  testthat::expect_s3_class(p, "ggplot")
})

# 5. fix_sample_alignment -------------------------------------------------
testthat::test_that("fix_sample_alignment aligns matrix and metadata by sample name", {
  d <- make_test_data()
  aligned <- fix_sample_alignment(t(d$numeric_data), d$metadata)
  testthat::expect_equal(colnames(aligned$data_matrix), rownames(aligned$metadata))
})

testthat::test_that("fix_sample_alignment rejects non-data.frame metadata", {
  d <- make_test_data()
  testthat::expect_error(fix_sample_alignment(t(d$numeric_data), list(a = 1)), "data.frame")
})

# 6. create_volcano_plot_labeled ------------------------------------------
testthat::test_that("create_volcano_plot_labeled works with and without classification", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  cls <- classify_lipids(colnames(norm))
  diff <- perform_differential_analysis(norm, d$metadata, "Sample Group", method = "limma")
  res <- diff$results[[1]]

  p1 <- create_volcano_plot_labeled(res)
  testthat::expect_s3_class(p1, "ggplot")

  p2 <- create_volcano_plot_labeled(res, classification_data = cls, color_by = "LipidGroup")
  testthat::expect_s3_class(p2, "ggplot")
})

testthat::test_that("create_volcano_plot_labeled requires logFC and adj.P.Val columns", {
  testthat::expect_error(
    create_volcano_plot_labeled(data.frame(x = 1:3)),
    "logFC"
  )
})

# 7. create_lipid_expression_barplot --------------------------------------
testthat::test_that("create_lipid_expression_barplot returns single plot or list", {
  d <- make_test_data()
  lipid <- colnames(d$numeric_data)[1]

  p_single <- create_lipid_expression_barplot(
    d$numeric_data, d$metadata,
    selected_lipids = lipid, group_column = "Sample Group"
  )
  testthat::expect_s3_class(p_single, "ggplot")

  p_multi <- create_lipid_expression_barplot(
    d$numeric_data, d$metadata,
    selected_lipids = colnames(d$numeric_data)[1:3], group_column = "Sample Group"
  )
  testthat::expect_type(p_multi, "list")
  testthat::expect_length(p_multi, 3)
})

testthat::test_that("create_lipid_expression_barplot rejects empty selected_lipids", {
  d <- make_test_data()
  testthat::expect_error(
    create_lipid_expression_barplot(d$numeric_data, d$metadata, selected_lipids = character(0)),
    "selected_lipids"
  )
})

testthat::test_that("create_lipid_expression_barplot errors when no lipids match", {
  d <- make_test_data()
  testthat::expect_error(
    create_lipid_expression_barplot(d$numeric_data, d$metadata, selected_lipids = "Not_A_Real_Lipid"),
    "None of the selected lipids"
  )
})

# 8. create_pca_plot_with_ellipses / create_plsda_plot_with_ellipses ------
testthat::test_that("create_pca_plot_with_ellipses works for all ellipse types", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  pca <- perform_pca(norm, d$metadata, "Sample Group")

  for (et in c("none", "confidence", "visual")) {
    p <- create_pca_plot_with_ellipses(pca$pca_data, pca$variance_explained, ellipse_type = et)
    testthat::expect_s3_class(p, "ggplot")
  }
})

testthat::test_that("create_pca_plot_with_ellipses validates pca_data columns", {
  testthat::expect_error(
    create_pca_plot_with_ellipses(data.frame(x = 1), c(50, 30)),
    "PC1"
  )
})

testthat::test_that("create_plsda_plot_with_ellipses works for all ellipse types", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  plsda <- perform_plsda(norm, d$metadata, "Sample Group")

  for (et in c("none", "confidence", "visual")) {
    p <- create_plsda_plot_with_ellipses(plsda$scores_data, ellipse_type = et)
    testthat::expect_s3_class(p, "ggplot")
  }
})

testthat::test_that("create_plsda_plot_with_ellipses validates plsda_data columns", {
  testthat::expect_error(
    create_plsda_plot_with_ellipses(data.frame(x = 1)),
    "Comp1"
  )
})

# 9. create_enrichment_dotplot / create_enrichment_barplot ----------------
testthat::test_that("enrichment plots work and handle empty results", {
  d <- make_test_data()
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  cls <- classify_lipids(colnames(norm))
  diff <- perform_differential_analysis(norm, d$metadata, "Sample Group", method = "limma")
  enrich <- perform_enrichment_analysis(diff$results, cls, min_set_size = 3)
  grp_enrich <- enrich[[1]][["LipidGroup"]]

  if (!is.null(grp_enrich) && nrow(grp_enrich) > 0) {
    p_dot <- create_enrichment_dotplot(grp_enrich)
    testthat::expect_s3_class(p_dot, "ggplot")
    p_bar <- create_enrichment_barplot(grp_enrich)
    testthat::expect_s3_class(p_bar, "ggplot")
  }

  p_empty <- create_enrichment_dotplot(data.frame())
  testthat::expect_s3_class(p_empty, "ggplot")
})

testthat::test_that("enrichment plots reject non-data.frame input", {
  testthat::expect_error(create_enrichment_dotplot(list(a = 1)), "data.frame")
  testthat::expect_error(create_enrichment_barplot(list(a = 1)), "data.frame")
})
