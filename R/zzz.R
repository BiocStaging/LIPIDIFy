# Register global variables to suppress R CMD check NOTEs from NSE usage
# in ggplot2, dplyr, tidyr, and related packages.
utils::globalVariables(c(
  ".data",
  "Abundance", "ColorGroup", "Comp1", "Comp2", "Direction",
  "Fill", "Group", "Intensity", "Label", "Lipid",
  "NES", "PC1", "PC2", "Sample", "Significance",
  "logFC", "neg_log10_pval", "padj", "pathway", "size"
))
