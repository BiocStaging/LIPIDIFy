# =============================================================================
# tests/testthat/test-normalization.R
#
# IMPORTANT: this file belongs in tests/testthat/, NOT in R/.
# Placing it in R/ causes devtools::load_all() to source it and fail because
# testthat::test_that() at the top level is run immediately during sourcing.
#
# Run with:
#   devtools::test()
# =============================================================================

make_matrix <- function(n_samples = 6, n_lipids = 20, seed = 42) {
  set.seed(seed)
  m <- matrix(stats::rlnorm(n_samples * n_lipids, meanlog = 8, sdlog = 1),
    nrow = n_samples, ncol = n_lipids
  )
  rownames(m) <- paste0("Sample_", seq_len(n_samples))
  colnames(m) <- paste0("Lipid_", seq_len(n_lipids))
  m
}

# 1. get_normalization_methods --------------------------------------------
testthat::test_that("get_normalization_methods returns character vector", {
  m <- get_normalization_methods()
  testthat::expect_type(m, "character")
  testthat::expect_true(length(m) > 0)
})

testthat::test_that("expected methods present", {
  m <- get_normalization_methods()
  for (x in c(
    "TIC", "PQN", "Quantile", "VSN", "Log2Median", "Median", "Mean",
    "Log2", "Log10", "Sqrt", "None"
  )) {
    testthat::expect_true(x %in% m, info = paste("Missing:", x))
  }
})

testthat::test_that("Min method has been removed (Mean is the correct scale method)", {
  testthat::expect_false("Min" %in% get_normalization_methods())
})

# 2. normalize_tic --------------------------------------------------------
testthat::test_that("normalize_tic: same dims", {
  m <- make_matrix()
  r <- normalize_tic(m)
  testthat::expect_equal(dim(r), dim(m))
})
testthat::test_that("normalize_tic: row sums equal", {
  m <- make_matrix()
  r <- normalize_tic(m)
  testthat::expect_equal(length(unique(round(rowSums(r), 4))), 1L)
})

# 3. normalize_pqn --------------------------------------------------------
testthat::test_that("normalize_pqn: same dims + rownames", {
  m <- make_matrix()
  r <- normalize_pqn(m)
  testthat::expect_equal(dim(r), dim(m))
  testthat::expect_equal(rownames(r), rownames(m))
})

# 4. normalize_quantile ---------------------------------------------------
testthat::test_that("normalize_quantile: same dims", {
  m <- make_matrix()
  r <- normalize_quantile(m)
  testthat::expect_equal(dim(r), dim(m))
})
testthat::test_that("normalize_quantile: identical sorted distributions (expected behaviour)", {
  m <- make_matrix()
  r <- normalize_quantile(m)
  sc <- apply(r, 1, sort)
  for (i in seq(2, ncol(sc))) {
    testthat::expect_equal(sc[, 1], sc[, i],
      tolerance = 1e-8,
      info = paste("Sample", i, "differs from sample 1")
    )
  }
})
testthat::test_that("normalize_quantile: warns for 1 sample", {
  testthat::expect_warning(normalize_quantile(make_matrix(n_samples = 1)))
})

# 5. normalize_log2median -------------------------------------------------
testthat::test_that("normalize_log2median: same dims", {
  m <- make_matrix()
  r <- normalize_log2median(m)
  testthat::expect_equal(dim(r), dim(m))
})
# 6. normalize_median -----------------------------------------------------
testthat::test_that("normalize_median: same dims", {
  m <- make_matrix()
  r <- normalize_median(m)
  testthat::expect_equal(dim(r), dim(m))
})
testthat::test_that("normalize_median: all sample medians equal", {
  m <- make_matrix()
  r <- normalize_median(m)
  testthat::expect_equal(length(unique(round(apply(r, 1, stats::median), 4))), 1L)
})

# 7. normalize_mean -------------------------------------------------------
testthat::test_that("normalize_mean: same dims", {
  m <- make_matrix()
  r <- normalize_mean(m)
  testthat::expect_equal(dim(r), dim(m))
})
testthat::test_that("normalize_mean: all sample means equal", {
  m <- make_matrix()
  r <- normalize_mean(m)
  testthat::expect_equal(length(unique(round(rowMeans(r), 4))), 1L)
})
testthat::test_that("normalize_median vs normalize_mean differ on skewed data", {
  # On symmetric (log-normal) data the two methods look similar — that is
  # expected and NOT a bug. On strongly skewed data they must differ.
  set.seed(1)
  sk <- matrix(c(1, 1, 1, 1, 1, 500, 2, 2, 2, 2, 2, 600), nrow = 2, byrow = TRUE)
  testthat::expect_false(
    isTRUE(all.equal(normalize_median(sk),
      normalize_mean(sk),
      tolerance = 1e-6
    )),
    info = "Median and Mean normalizations must differ on skewed data"
  )
})

# 8. apply_normalizations -------------------------------------------------
testthat::test_that("apply_normalizations None is identity", {
  m <- make_matrix()
  testthat::expect_equal(apply_normalizations(m, "None"), m)
})
testthat::test_that("apply_normalizations chains correctly", {
  m <- make_matrix()
  testthat::expect_equal(
    apply_normalizations(m, c("TIC", "Log2")),
    log2(normalize_tic(m) + 1)
  )
})
testthat::test_that("apply_normalizations preserves rownames", {
  m <- make_matrix()
  testthat::expect_equal(rownames(apply_normalizations(m, "TIC")), rownames(m))
})
testthat::test_that("apply_normalizations passes through unknown method", {
  m <- make_matrix()
  testthat::expect_equal(apply_normalizations(m, "XYZ_UNKNOWN"), m)
})

# 9. classify_lipids / determine_saturation --------------------------------
testthat::test_that("classify_lipids returns correct columns", {
  r <- classify_lipids(c("PC 16:0_18:1", "PE 18:0_20:4"))
  testthat::expect_true(all(c("Lipid", "LipidGroup", "LipidType", "Saturation") %in% colnames(r)))
  testthat::expect_equal(nrow(r), 2L)
})
testthat::test_that("classify_lipids assigns groups correctly", {
  r <- classify_lipids(c("PC 16:0_18:1", "TG 16:0_18:1_20:4", "SM 18:1"))
  testthat::expect_equal(r$LipidGroup[r$Lipid == "PC 16:0_18:1"], "Glycerophospholipids")
  testthat::expect_equal(r$LipidGroup[r$Lipid == "TG 16:0_18:1_20:4"], "Glycerolipids")
  testthat::expect_equal(r$LipidGroup[r$Lipid == "SM 18:1"], "Sphingolipids")
})
testthat::test_that("determine_saturation is correct", {
  testthat::expect_equal(determine_saturation("PC 18:0_18:0"), "SFA")
  testthat::expect_equal(determine_saturation("PC 16:0_18:1"), "MUFA")
  testthat::expect_equal(determine_saturation("LPC 20:4"), "PUFA")
  testthat::expect_equal(determine_saturation("TG 16:0_18:1_20:4"), "PUFA")
  testthat::expect_equal(determine_saturation("PE(P-18:0_20:5)"), "PUFA")
})

# 10. load_lipidomics_data_from_df ----------------------------------------
testthat::test_that("load_lipidomics_data_from_df basic structure", {
  res <- load_lipidomics_data_from_df(generate_example_data())
  testthat::expect_true(is.matrix(res$numeric_data))
  testthat::expect_true(is.data.frame(res$metadata))
  testthat::expect_equal(nrow(res$numeric_data), nrow(res$metadata))
})
testthat::test_that("load_lipidomics_data_from_df removes PBQC", {
  ex <- generate_example_data()
  ex$`Sample Group`[1] <- "PBQC"
  res <- load_lipidomics_data_from_df(ex)
  testthat::expect_false(any(grepl("PBQC", res$metadata$`Sample Group`, ignore.case = TRUE)))
})

# 11. create_default_contrasts --------------------------------------------
testthat::test_that("create_default_contrasts: all pairwise for 3 groups", {
  testthat::expect_equal(length(create_default_contrasts(c("A", "B", "C"))), 3L)
})
testthat::test_that("create_default_contrasts: empty for 1 group", {
  testthat::expect_length(create_default_contrasts("OnlyGroup"), 0L)
})
testthat::test_that("create_default_contrasts: sanitizes spaces in names", {
  ctrs <- create_default_contrasts(c("Group A", "Group B"))
  testthat::expect_length(ctrs, 1L)
  # The contrast uses limma format "TokenB - TokenA"; the ' - ' separator is
  # intentional and always present. Check that the group-name *tokens* have
  # had their spaces removed by make.names(), not the separator itself.
  tokens <- trimws(unlist(strsplit(ctrs, " - ")))
  testthat::expect_false(any(grepl(" ", tokens)),
    info = "Group-name tokens must not contain spaces after make.names() sanitization"
  )
})
