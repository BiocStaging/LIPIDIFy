# =============================================================================
# tests/testthat/test-utils.R
# Coverage for R/utils.R not already exercised by test-normalization.R:
# load_lipidomics_data, imputation, batch correction, classification I/O,
# matrix orientation, and shared validators.
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

# 1. load_lipidomics_data ------------------------------------------------
testthat::test_that("load_lipidomics_data reads a CSV and separates metadata/numeric", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(
    data.frame(
      "Sample Name" = c("S1", "S2"), "Sample Group" = c("A", "B"),
      "PC 16:0" = c(1000, 2000), check.names = FALSE
    ),
    tmp,
    row.names = FALSE
  )
  res <- load_lipidomics_data(tmp)
  testthat::expect_equal(dim(res$numeric_data), c(2L, 1L))
  testthat::expect_equal(nrow(res$metadata), 2L)
})

testthat::test_that("load_lipidomics_data errors on missing file", {
  testthat::expect_error(load_lipidomics_data("does/not/exist.csv"), "not found")
})

testthat::test_that("load_lipidomics_data errors on non-character file_path", {
  testthat::expect_error(load_lipidomics_data(123), "file_path")
})

testthat::test_that("load_lipidomics_data removes PBQC rows", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(
    data.frame(
      "Sample Name" = c("S1", "S2", "S3"),
      "Sample Group" = c("A", "B", "PBQC"),
      "PC 16:0" = c(1000, 2000, 3000), check.names = FALSE
    ),
    tmp,
    row.names = FALSE
  )
  res <- load_lipidomics_data(tmp)
  testthat::expect_equal(nrow(res$metadata), 2L)
})

# 2. .validate_metadata / .validate_design (internal) ---------------------
testthat::test_that(".validate_metadata errors on non-data.frame or missing column", {
  testthat::expect_error(LIPIDIFy:::.validate_metadata(list(a = 1)), "data.frame")
  testthat::expect_error(
    LIPIDIFy:::.validate_metadata(data.frame(a = 1), "missing_col"),
    "missing_col"
  )
  testthat::expect_true(LIPIDIFy:::.validate_metadata(data.frame(a = 1)))
})

testthat::test_that(".validate_design aligns by rownames and errors on mismatch", {
  m <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("S1", "S2", "S3")))
  design <- matrix(1, nrow = 3, ncol = 1, dimnames = list(c("S3", "S1", "S2"), "grp"))
  aligned <- LIPIDIFy:::.validate_design(design, m)
  testthat::expect_equal(rownames(aligned), c("S1", "S2", "S3"))

  bad_design <- matrix(1, nrow = 2, ncol = 1, dimnames = list(c("Sx", "Sy"), "grp"))
  testthat::expect_error(LIPIDIFy:::.validate_design(bad_design, m), "design")
})

# 4. classification I/O ----------------------------------------------------
testthat::test_that("load_custom_classification / export_classification round-trip", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(
    data.frame(Lipid = c("PC 16:0_18:1", "PE 18:0"), Class = c("A", "B")),
    tmp,
    row.names = FALSE
  )
  cls <- load_custom_classification(tmp)
  testthat::expect_equal(colnames(cls)[1], "Lipid")

  tmp2 <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp2), add = TRUE)
  export_classification(cls, tmp2)
  testthat::expect_true(file.exists(tmp2))
})

testthat::test_that("load_custom_classification requires a Lipid column", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(NotLipid = "x"), tmp, row.names = FALSE)
  testthat::expect_error(load_custom_classification(tmp), "Lipid")
})

testthat::test_that("load_custom_classification errors on missing file", {
  testthat::expect_error(load_custom_classification("nope.csv"), "not found")
})

testthat::test_that("export_classification requires a data.frame", {
  testthat::expect_error(export_classification(list(a = 1), tempfile()), "data.frame")
})

testthat::test_that("load_custom_enrichment_sets parses Lipid/Set_Name pairs", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(
    data.frame(
      Lipid = c("PC 16:0_18:1", "PE 18:0"),
      Set_Name = c("SetA", "SetA")
    ),
    tmp,
    row.names = FALSE
  )
  sets <- load_custom_enrichment_sets(tmp)
  testthat::expect_named(sets, "SetA")
  testthat::expect_length(sets$SetA, 2)
})

testthat::test_that("load_custom_enrichment_sets requires Lipid and Set_Name columns", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(data.frame(Lipid = "x"), tmp, row.names = FALSE)
  testthat::expect_error(load_custom_enrichment_sets(tmp), "Set_Name")
})

# 5. imputation -------------------------------------------------------------
testthat::test_that("impute_missing_values half_min/min/zero/mean/median remove all NAs", {
  m <- make_matrix()
  m[1, 1] <- NA
  m[2, 3] <- NA
  for (method in c("half_min", "min", "zero", "mean", "median")) {
    imputed <- impute_missing_values(m, method = method)
    testthat::expect_equal(sum(is.na(imputed)), 0L, info = method)
    testthat::expect_equal(dim(imputed), dim(m))
  }
})

testthat::test_that("impute_missing_values returns data unchanged when nothing is missing", {
  m <- make_matrix()
  testthat::expect_message(r <- impute_missing_values(m), "No missing values")
  testthat::expect_equal(r, m)
})

testthat::test_that("impute_missing_values rejects invalid data_matrix", {
  testthat::expect_error(impute_missing_values("not a matrix"), "data_matrix")
})

testthat::test_that("impute_missing_values knn seed does not leak into global RNG state", {
  testthat::skip_if_not_installed("impute")
  m <- make_matrix()
  m[1, 1] <- NA

  set.seed(123)
  before <- runif(1)

  set.seed(123)
  invisible(impute_missing_values(m, method = "knn", seed = 999))
  after <- runif(1)

  testthat::expect_equal(before, after)
})

# 6. batch effect correction -------------------------------------------------
testthat::test_that("correct_batch_effects (limma) preserves dimensions and reports method_used", {
  d <- load_lipidomics_data_from_df(generate_example_data())
  md <- d$metadata
  md$Batch <- rep(c("B1", "B2"), times = 10)
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))

  corrected <- suppressWarnings(
    correct_batch_effects(norm, md, batch_column = "Batch", method = "limma")
  )
  testthat::expect_equal(dim(corrected), dim(norm))
  testthat::expect_equal(attr(corrected, "method_used"), "limma")
})

testthat::test_that("correct_batch_effects falls back from combat to limma when sva is unavailable", {
  testthat::skip_if(requireNamespace("sva", quietly = TRUE), "sva is installed; fallback path not triggered")
  d <- load_lipidomics_data_from_df(generate_example_data())
  md <- d$metadata
  md$Batch <- rep(c("B1", "B2"), times = 10)
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))

  corrected <- suppressWarnings(
    correct_batch_effects(norm, md, batch_column = "Batch", method = "combat")
  )
  testthat::expect_equal(attr(corrected, "method_used"), "limma")
})

testthat::test_that("correct_batch_effects errors on missing batch_column", {
  d <- load_lipidomics_data_from_df(generate_example_data())
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(
    correct_batch_effects(norm, d$metadata, batch_column = "Nope"),
    "batch_column"
  )
})

testthat::test_that("correct_batch_effects rejects non-data.frame metadata", {
  d <- load_lipidomics_data_from_df(generate_example_data())
  norm <- apply_normalizations(d$numeric_data, c("TIC", "Log2"))
  testthat::expect_error(
    correct_batch_effects(norm, list(a = 1), batch_column = "Batch"),
    "data.frame"
  )
})

# 7. apply_normalizations validity checks -----------------------------------
testthat::test_that("apply_normalizations rejects invalid data or methods", {
  m <- make_matrix()
  testthat::expect_error(apply_normalizations("not data", "TIC"), "data")
  testthat::expect_error(apply_normalizations(m, character(0)), "methods")
})

testthat::test_that("individual normalize_* functions reject non-numeric input", {
  bad <- data.frame(a = c("x", "y"), b = c("z", "w"))
  for (fn in list(
    normalize_tic, normalize_pqn, normalize_quantile,
    normalize_log2median, normalize_median, normalize_mean
  )) {
    testthat::expect_error(fn(bad), "numeric")
  }
  testthat::expect_error(normalize_tic(list(a = 1)), "matrix or data.frame")
})

testthat::test_that(".validate_numeric_matrix accepts numeric matrices and data frames", {
  m <- make_matrix()
  df <- as.data.frame(m)
  testthat::expect_true(LIPIDIFy:::.validate_numeric_matrix(m))
  testthat::expect_true(LIPIDIFy:::.validate_numeric_matrix(df))
})

# 8. get_imputation_descriptions / get_normalization_descriptions -----------
testthat::test_that("descriptions cover every method name", {
  testthat::expect_setequal(names(get_imputation_descriptions()), get_imputation_methods())
  testthat::expect_setequal(names(get_normalization_descriptions()), get_normalization_methods())
})
