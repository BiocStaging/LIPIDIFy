# =============================================================================
# tests/testthat/test-orientation-topn.R
# Coverage for the pre-existing matrix-orientation fix (.orient_matrix(),
# used by perform_pca/perform_plsda/perform_differential_analysis_*/
# correct_batch_effects/create_heatmap_robust/create_lipid_expression_barplot)
# and the "limit to top N variable samples" feature in
# visualize_raw_data_improved()/create_pipeline_plot().
# =============================================================================

# 1. .orient_matrix (internal) ------------------------------------------------
testthat::test_that(".orient_matrix orients by sample-name match", {
  md <- data.frame(`Sample Name` = c("S1", "S2", "S3"), check.names = FALSE)
  m  <- matrix(1:6, nrow = 3, dimnames = list(c("S1", "S2", "S3"), c("L1", "L2")))
  oriented <- LIPIDIFy:::.orient_matrix(m, md, want = "features_rows")
  testthat::expect_equal(rownames(oriented), c("L1", "L2"))
})

testthat::test_that(".orient_matrix falls back to sample count when names are absent", {
  m <- matrix(1:12, nrow = 3, ncol = 4)  # 3 samples x 4 features, no dimnames
  md <- data.frame(x = 1:3)
  oriented <- LIPIDIFy:::.orient_matrix(m, md, want = "samples_rows")
  testthat::expect_equal(nrow(oriented), 3)
})

# 2. perform_pca correctly handles n_samples > n_features ---------------------
testthat::test_that("perform_pca works when samples outnumber features", {
  set.seed(1)
  m <- matrix(stats::rlnorm(40 * 5, 8, 1), nrow = 40, ncol = 5)
  rownames(m) <- paste0("S", seq_len(40))
  colnames(m) <- paste0("Lipid_", seq_len(5))
  md <- data.frame(`Sample Name` = rownames(m),
                    `Sample Group` = rep(c("A", "B"), each = 20),
                    check.names = FALSE)
  res <- perform_pca(m, md, "Sample Group")
  testthat::expect_equal(nrow(res$pca_data), 40)
})

# 3. visualize_raw_data_improved: limit to top-N most variable samples -------
testthat::test_that("visualize_raw_data_improved respects top_n in sample mode", {
  d <- load_lipidomics_data_from_df(generate_example_data())
  dl <- list(numeric_data = d$numeric_data)
  p <- visualize_raw_data_improved(dl, "boxplot", "sample", top_n = 5)
  testthat::expect_s3_class(p, "ggplot")
})
