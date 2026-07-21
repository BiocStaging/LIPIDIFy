# Prepare the real-data example bundled at inst/extdata/ST001359_lipidomics.csv
#
# Source: NIH Common Fund Metabolomics Workbench, Study ST001359
#   "Monophasic lipidomics extraction in cancer cell lines"
#   Rodriguez Blanco G., Beatson Institute for Cancer Research
#   Project ID: PR000929, Analysis ID: AN002263
#   https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST001359
#   License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
#
# Design: HepG2 hepatocellular carcinoma cells, LC-MS (reversed phase,
# positive ion mode), 6 samples = 3 untreated controls vs. 3 cells treated
# with an SDC (syndecan) inhibitor. Reported values are the study's own
# "Peak area normalized" units (already normalized and log2-scaled by the
# submitters) -- they are used here unchanged, i.e. not renormalized or
# back-transformed, to preserve fidelity to the deposited data.
#
# This script is not run as part of the package build; it documents how
# inst/extdata/ST001359_lipidomics.csv was derived from the public API and
# is kept for provenance/reproducibility only.

raw <- jsonlite::fromJSON(
  "https://www.metabolomicsworkbench.org/rest/study/study_id/ST001359/data",
  simplifyVector = FALSE
)

lipid_names <- vapply(raw, function(x) x$metabolite_name, character(1))
sample_ids <- names(raw[[1]]$DATA)

mat <- vapply(raw, function(x) as.numeric(unlist(x$DATA)), numeric(length(sample_ids)))
rownames(mat) <- sample_ids
colnames(mat) <- lipid_names

sample_group <- ifelse(grepl("_C[0-9]+$", sample_ids), "Control", "SDC_i")

out <- data.frame(
  "Sample Name" = sample_ids,
  "Sample Group" = sample_group,
  mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(out, "inst/extdata/ST001359_lipidomics.csv", row.names = FALSE)
