# Register global variables to suppress R CMD check NOTEs from NSE usage
# in ggplot2, dplyr, tidyr, and related packages.
utils::globalVariables(c(
  ".data",
  "Abundance", "ColorGroup", "Comp1", "Comp2", "Direction",
  "Fill", "Group", "Intensity", "Label", "Lipid",
  "NES", "PC1", "PC2", "Sample", "Significance",
  "logFC", "neg_log10_pval", "padj", "pathway", "size"
))

# ---------------------------------------------------------------------------
# Package-level imports
# These directives are picked up by roxygen2 (devtools::document()) and
# written to NAMESPACE as importFrom() entries, eliminating R CMD check
# NOTEs about undefined global functions.
# ---------------------------------------------------------------------------

#' @importFrom utils head write.csv capture.output
#' @importFrom stats model.matrix median setNames
NULL
