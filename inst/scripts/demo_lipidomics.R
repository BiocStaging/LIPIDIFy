# End-to-end LIPIDIFy demo on a real public lipidomics dataset.
#
# Data: Metabolomics Workbench Study ST001359, "Monophasic lipidomics
# extraction in cancer cell lines" (Rodriguez Blanco G., Beatson Institute
# for Cancer Research; CC BY 4.0). HepG2 cells, Control vs. SDC-inhibitor
# treatment, 3 replicates each. See inst/extdata/ST001359_lipidomics.csv
# and data-raw/prepare_ST001359_lipidomics.R for provenance details.
#
# Run with: Rscript system.file("scripts", "demo_lipidomics.R", package = "LIPIDIFy")

library(LIPIDIFy)

file_path <- system.file("extdata", "ST001359_lipidomics.csv", package = "LIPIDIFy")

loaded <- load_lipidomics_data(
  file_path,
  metadata_columns = c("Sample Name", "Sample Group")
)

str(loaded$metadata)
dim(loaded$numeric_data)

# Classify lipids
lipid_names <- colnames(loaded$numeric_data)
classification <- classify_lipids(lipid_names)

head(classification)
table(classification$LipidGroup)

# Normalize (values are already log2-scaled by the source; median-centre only)
norm_mat <- apply_normalizations(loaded$numeric_data, c("Median"))

# PCA
pca_res <- perform_pca(norm_mat, loaded$metadata, group_column = "Sample Group")

p_pca <- create_pca_plot_with_ellipses(
  pca_data = pca_res$pca_data,
  variance_explained = pca_res$variance_explained,
  ellipse_type = "confidence"
)
print(p_pca)

dir.create("demo_output", showWarnings = FALSE)
ggplot2::ggsave("demo_output/PCA_normalized.png", p_pca, width = 8, height = 6, dpi = 300)

# Differential analysis: Control vs SDC_i
contrasts <- create_default_contrasts(unique(loaded$metadata$`Sample Group`))

diff_res <- perform_differential_analysis(
  data_matrix = norm_mat,
  metadata = loaded$metadata,
  group_column = "Sample Group",
  contrasts_list = contrasts,
  method = "limma"
)

first_contrast_name <- names(diff_res$results)[1]
de_table <- diff_res$results[[first_contrast_name]]
head(de_table)
summary(de_table$adj.P.Val < 0.05)

# Volcano plot
p_volcano <- create_volcano_plot_labeled(
  results = de_table,
  classification_data = classification,
  color_by = "LipidGroup"
)
print(p_volcano)
ggplot2::ggsave("demo_output/volcano.png", p_volcano, width = 8, height = 6, dpi = 300)

# Enrichment analysis
enrich <- perform_enrichment_analysis(
  results_list = diff_res$results,
  classification_data = classification,
  min_set_size = 5,
  max_set_size = 500,
  custom_sets = NULL
)

first_contrast_enrich <- enrich[[first_contrast_name]]

# Pick one enrichment type with non-empty results
enrich_df <- NULL
for (et in names(first_contrast_enrich)) {
  if (nrow(first_contrast_enrich[[et]]) > 0) {
    enrich_df <- first_contrast_enrich[[et]]
    enrich_type <- et
    break
  }
}

if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
  p_enrich <- create_enrichment_dotplot(
    enrichment_data = enrich_df,
    title = paste("Enrichment:", first_contrast_name, "-", enrich_type),
    max_pathways = 15
  )
  print(p_enrich)
  ggplot2::ggsave(
    "demo_output/enrichment_first_contrast.png",
    p_enrich,
    width = 8, height = 6, dpi = 300
  )
} else {
  message("No non-empty enrichment results to plot.")
}
