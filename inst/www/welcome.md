# Welcome to LIPIDIFy

## Overview

**LIPIDIFy** is a comprehensive R package for lipidomics data analysis, providing both a user-friendly Shiny interface and programmatic access through R functions. This package is designed to facilitate:

- Data upload and preprocessing
- Multiple normalization strategies
- Differential abundance analysis (limma and EdgeR)
- Enrichment analysis
- Data visualization (PCA, PLS-DA, heatmaps, volcano plots)
- Lipid classification and annotation

## Getting Started

### For Biologists (Shiny Interface)

Navigate through the tabs on the left to:
1. **Data Upload**: Load your CSV data or use example data
2. **Lipid Classification**: View/download automatic classification or upload custom classification
3. **Raw Data Visualization**: Explore your raw data
4. **Normalization**: Compare and apply normalization methods
5. **Differential Analysis**: Identify differentially abundant lipids
6. **Enrichment Analysis**: Find enriched lipid classes

### For Bioinformaticians (R Functions)

All functions are available for programmatic use:

```r
# Load data
data <- load_lipidomics_data("your_data.csv")

# Normalize
normalized <- apply_normalizations(data$numeric_data, methods = c("TIC", "Log2"))

# Differential analysis
results <- perform_differential_analysis(normalized, data$metadata)

# Classification
classification <- classify_lipids(colnames(data$numeric_data))
```

---

## Input Data Structure

Your input CSV file should have the following structure:

### Example Data Format

| Sample Name | Sample Group | Tumour ID | Weight (mg) | PC 16:0_18:1 | PE 18:0_20:4 | TG 16:0_18:1_20:4 | ... |
|-------------|--------------|-----------|-------------|--------------|--------------|-------------------|-----|
| Sample_1    | GroupA       | 1         | 45.2        | 125000       | 89000        | 234000            | ... |
| Sample_2    | GroupA       | 1         | 48.7        | 118000       | 92000        | 221000            | ... |
| Sample_3    | GroupB       | 2         | 52.1        | 135000       | 95000        | 245000            | ... |
| Sample_4    | GroupB       | 2         | 46.8        | 142000       | 88000        | 251000            | ... |

**Requirements:**
- First columns: Metadata (Sample Name, Sample Group are required)
- Remaining columns: Lipid abundances (numeric values)
- Lipid names as column headers (e.g., "PC 16:0_18:1")
- CSV format with comma separation

---

## Lipid Classification Structure

The package automatically classifies lipids based on their names. You can also provide custom classifications.

### Automatic Classification Example

When you load data, lipids are automatically classified:

| Lipid | LipidGroup | LipidType | Saturation |
|-------|------------|-----------|------------|
| PC 16:0_18:1 | Glycerophospholipids | Phosphatidylcholine | MUFA |
| PE 18:0_20:4 | Glycerophospholipids | Phosphatidylethanolamine | PUFA |
| TG 16:0_18:1_20:4 | Glycerolipids | Triacylglycerol | PUFA |
| SM 16:0 | Sphingolipids | Sphingomyelin | SFA |
| Cer 18:1 | Sphingolipids | Ceramide | MUFA |

**Classification Levels:**
- **LipidGroup**: Main lipid category (e.g., Glycerophospholipids, Sphingolipids)
- **LipidType**: Specific lipid class (e.g., Phosphatidylcholine, Ceramide)
- **Saturation**: Fatty acid saturation (SFA, MUFA, PUFA, or Unclassified)

### Custom Classification Format

You can provide your own classification with flexible hierarchical levels:

**Minimum format (Lipid + 1 classification level):**

| Lipid | Class1 |
|-------|--------|
| PC 16:0_18:1 | Phospholipids |
| PE 18:0_20:4 | Phospholipids |
| TG 16:0_18:1_20:4 | Neutral Lipids |

**Extended format (Lipid + multiple classification levels):**

| Lipid | MainClass | SubClass | FunctionalGroup | Saturation |
|-------|-----------|----------|-----------------|------------|
| PC 16:0_18:1 | Phospholipids | Glycerophospholipids | Membrane Components | MUFA |
| PE 18:0_20:4 | Phospholipids | Glycerophospholipids | Membrane Components | PUFA |
| TG 16:0_18:1_20:4 | Neutral Lipids | Glycerolipids | Energy Storage | PUFA |

**Requirements for custom classification:**
- First column must be named "Lipid"
- Subsequent columns can have any names (Class1, Class2, MainClass, etc.)
- Lipid names must match those in your data file
- CSV format

Upload custom classifications in the **Lipid Classification** tab.

---

## Analysis Workflow

```
Data Upload → Classification → Raw Visualization → Normalization →
→ Normalized Visualization → Differential Analysis → Results Visualization → Enrichment
```

### Key Features

#### Normalization Methods
- TIC (Total Ion Current)
- PQN (Probabilistic Quotient Normalization)
- Quantile, VSN, Median, Mean
- Log transformations (Log2, Log10)
- Pipeline comparison tools

#### Differential Analysis
- **limma**: Linear models for microarray/omics data
- **EdgeR**: Robust method originally for RNA-seq, adaptable to lipidomics
- Flexible contrast selection
- Multiple testing correction (FDR)

#### Visualization
- PCA and PLS-DA with confidence ellipses
- Heatmaps with clustering
- Volcano plots
- Barplots for individual lipid expression
- Enrichment dotplots and barplots

---

## Citation

If you use LIPIDIFy in your research, please cite:

*[Citation information will be added upon publication]*

---

## Support

For questions, issues, or feature requests, please visit:
- GitHub: [Repository link]
- Documentation: Run `?function_name` in R for help on specific functions

---

**Ready to start?** Click on **Data Upload** in the sidebar to begin your analysis!
