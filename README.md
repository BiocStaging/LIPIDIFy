# LIPIDIFy

**LIPIDIFy** is a comprehensive R package for lipidomics data analysis, designed to support both:

- **Experimental biologists**, via an interactive Shiny application
- **Bioinformaticians and data scientists**, via a fully scriptable R workflow

The package provides an end-to-end pipeline from raw lipidomics data to biological interpretation.

---

## Features

- Flexible CSV data import with automatic metadata detection
- **Automatic lipid classification** by class, subclass, and fatty-acid saturation (SFA/MUFA/PUFA)
- Custom classification upload (any number of classification columns) with reset to automatic
- Ten **normalization strategies**: TIC, PQN, Quantile, VSN, Median, Mean, Log2, Log10, Sqrt, None
- Pipeline builder: chain methods in any order and compare two pipelines side-by-side
- **Differential abundance analysis** (limma, EdgeR) with automatic pairwise contrasts
- **Enrichment analysis** of lipid classes or custom pathway sets (fgsea)
- Rich **visualizations**:
  - PCA and PLS-DA with confidence ellipses and sample labels
  - Heatmaps of significant features
  - Volcano plots coloured by lipid class
  - Per-sample and per-lipid barplots and boxplots
  - Density plots and violin plots
- Sample group filtering throughout the app
- Interactive **Shiny interface** (no coding required)
- Reproducible, modular R functions for scripted workflows
- Automated HTML or PDF report generation

---

## Installation

### Bioconductor (release)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LIPIDIFy")
library(LIPIDIFy)
```

### Development version

```r
BiocManager::install("fayrouzhammal/LIPIDIFy")
library(LIPIDIFy)
```

---

## Use Case 1 - Shiny App (no coding required)

Launch the interactive application:

```r
launch_lipidomics_app()
```

### Shiny workflow

1. **Data Upload** - upload a CSV file or load the built-in example dataset
2. **Lipid Classification** - automatic classification from lipid names, or upload a custom classification CSV; reset to automatic at any time
3. **Raw Data Visualization** - inspect intensity distributions before processing; filter by sample group; choose boxplot, violin, density, or histogram
4. **Normalization** - build and compare two normalization pipelines side-by-side; apply the best one
5. **Normalized Data Visualization** - verify normalization quality with boxplots, PCA, or PLS-DA; filter groups; add confidence ellipses
6. **Lipid Expression** - per-lipid abundance barplots across samples, sorted by group
7. **Differential Analysis** - limma or EdgeR; all pairwise contrasts generated automatically; FDR correction
8. **Results Visualization** - volcano plots (coloured by lipid class) or heatmaps
9. **Enrichment Analysis** - enrichment of lipid classes or custom pathway sets; dot plots or bar plots
10. **Generate Report** - export all figures and tables as HTML or PDF

All plots and tables can be downloaded directly from the interface.

---

## Use Case 2 - Programmatic R Workflow

### Load lipidomics data

```r
data <- load_lipidomics_data("your_data.csv")
# Returns a list with:
#   $numeric_data  - matrix of lipid abundances (samples x lipids)
#   $metadata      - data frame of sample metadata
#   $data          - original full data frame
```

### Classify lipids

```r
classification <- classify_lipids(colnames(data$numeric_data))
# Returns LipidGroup, LipidType, Saturation for each lipid
```

### Normalize data

```r
# View available methods
get_normalization_methods()
# [1] "TIC"      "PQN"      "Quantile" "VSN"      "Median"   "Mean"
# [7] "Log2"     "Log10"    "Sqrt"     "None"

# Apply a pipeline (methods chained in order)
normalized <- apply_normalizations(
  data$numeric_data,
  methods = c("TIC", "Log2")
)
```

### Differential analysis

```r
results <- perform_differential_analysis(
  data_matrix    = normalized,
  metadata       = data$metadata,
  group_column   = "Sample Group",   # column containing group labels
  contrasts_list = NULL,             # NULL = all pairwise contrasts
  method         = "limma"           # or "edger"
)
# Returns a list with $results (named list of data frames, one per contrast)
# Each data frame has columns: logFC, AveExpr, t, P.Value, adj.P.Val, B
```

### Enrichment analysis

```r
enrichment <- perform_enrichment_analysis(
  results_list       = results$results,   # named list from perform_differential_analysis
  classification_data = classification,
  min_set_size       = 5,
  max_set_size       = 500,
  custom_sets        = NULL              # optional: named list of lipid vectors
)
```

### Visualize results

```r
# Volcano plot coloured by lipid class
create_volcano_plot_labeled(
  results         = results$results[["GroupB - GroupA"]],
  classification_data = classification,
  color_by        = "LipidGroup"
)

# PCA plot
pca_res <- perform_pca(normalized, data$metadata, "Sample Group")
create_pca_plot_with_ellipses(
  pca_data           = pca_res$pca_data,
  variance_explained = pca_res$variance_explained,
  ellipse_type       = "confidence"   # "none", "confidence", or "visual"
)
```

---

## Input Data Format

### File 1 - Main lipidomics data (required)

Format: **CSV**. Rows = samples. First columns = metadata. Remaining columns = lipid abundances.

| Sample Name | Sample Group | Tumour ID | Weight (mg) | PC 16:0_18:1 | PE 18:0_20:4 | TG 16:0_18:1_20:4 |
|-------------|-------------|-----------|-------------|--------------|--------------|-------------------|
| Sample_1 | Control | T001 | 45.2 | 125000 | 89000 | 234000 |
| Sample_2 | Control | T002 | 48.7 | 118000 | 92000 | 221000 |
| Sample_3 | Treatment | T003 | 52.1 | 135000 | 95000 | 245000 |

**Requirements:**
- `Sample Name` and `Sample Group` are required metadata columns
- Any number of additional metadata columns are supported (Tumour ID, Weight, etc.)
- All lipid abundance values must be positive numbers
- Lipid names should follow standard notation (see below)
- PBQC/QC pool rows (Sample Group = "PBQC") are removed automatically

**Standard lipid name notation:**

| Example | Class |
|---------|-------|
| `PC 16:0_18:1` | Phosphatidylcholine |
| `PE 18:0_20:4` | Phosphatidylethanolamine |
| `TG 16:0_18:1_20:4` | Triacylglycerol |
| `LPC 18:2` | Lysophosphatidylcholine |
| `SM 18:1` | Sphingomyelin |
| `Cer 16:0` | Ceramide |
| `CE 18:1` | Cholesteryl ester |

### File 2 - Custom lipid classification (optional)

Format: **CSV**. First column must be named `Lipid`. Any number of additional classification columns.

| Lipid | LipidClass | BiologicalRole |
|-------|-----------|----------------|
| PC 16:0_18:1 | Phospholipid | Structural |
| TG 16:0_18:1_20:4 | Glycerolipid | Energy Storage |
| Cer 16:0 | Sphingolipid | Apoptosis |

> **Important:** Lipid names must match your data file **exactly** (case-sensitive,
> character-for-character). A mismatch will result in lipids appearing as "Unknown"
> in volcano plots. The app shows a diagnostic message with example names from both
> sides to help identify mismatches.

Each classification column name automatically becomes a coloring option in volcano plots
and an enrichment category.

### File 3 - Custom enrichment sets (optional)

Format: **CSV** with exactly two columns: `Lipid` and `Set_Name`.
A lipid can appear in multiple rows to belong to multiple sets.

| Lipid | Set_Name |
|-------|---------|
| PC 16:0_18:1 | Membrane_Phospholipids |
| PE 18:0_20:4 | Membrane_Phospholipids |
| PE 18:0_20:4 | AA_Derived_Lipids |
| TG 16:0_18:1_18:1 | Energy_Storage |
| Cer 16:0 | Apoptosis_Ceramides |

---

## Analysis Workflow

```
Data Upload
  -> Lipid Classification (automatic or custom)
  -> Raw Data Visualization (QC)
  -> Normalization (pipeline builder + comparison)
  -> Normalized Data Visualization (PCA, boxplots)
  -> Lipid Expression (per-lipid barplots)
  -> Differential Analysis (limma / EdgeR)
  -> Results Visualization (volcano, heatmap)
  -> Enrichment Analysis (lipid class or custom sets)
  -> Report Generation (HTML or PDF)
```

---

## Normalization Methods

| Method | Description | Common use |
|--------|-------------|-----------|
| TIC | Divides each sample by its total ion current, rescaled to global mean | Variable injection amounts |
| PQN | Probabilistic Quotient Normalization; robust to large fold-changes | General purpose |
| Quantile | Forces identical distributions across all samples | Removes all between-sample distribution differences |
| VSN | Variance Stabilizing (simplified): log2 + median centering | Mean-variance dependent data |
| Median | Scales each sample to the global median | Simple and robust |
| Mean | Scales each sample to the global mean | Similar to Median |
| Log2 | Log base-2 transformation | Almost always applied as last step |
| Log10 | Log base-10 transformation | Alternative to Log2 |
| Sqrt | Square-root transformation | Milder variance stabilization |
| None | No transformation | Data already normalized |

Methods are applied in the order selected. Recommended pipelines: **TIC + Log2** or **PQN + Log2**.

> **Note on Quantile normalization:** after applying Quantile normalization, all per-sample
> boxplots will look nearly identical. This is the **expected, correct** behaviour, not a bug.

---

## Documentation

```r
# Function-level help
?launch_lipidomics_app
?apply_normalizations
?perform_differential_analysis

# Full vignette
browseVignettes("LIPIDIFy")
```

---

## Citation

If you use **LIPIDIFy** in your research, please cite:

> Citation information will be added upon publication.

---

## Support

- Bug reports and feature requests: GitHub **Issues**
- Contributions welcome via pull requests
