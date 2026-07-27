# LIPIDIFy 0.99.3

## New features
- Added `normalize_vsn()`, a genuine variance-stabilising normalization that
  calls `vsn::justvsn()` from the Bioconductor **vsn** package (Huber et al.
  2002). `"VSN"` is now a separate accepted method in
  `get_normalization_methods()`, `apply_normalizations()` and the Shiny
  normalization pipeline builder.
- `vsn` added to `Suggests`; it is loaded dynamically via
  `requireNamespace()`, matching how the package already handles `impute`
  and `sva`.

## Behaviour changes
- **`"VSN"` and `"Log2Median"` are two distinct methods.** Previously the
  README, NEWS and Shiny help text described a `"VSN"` method that did not
  exist in the code: the implemented method was `"Log2Median"` (a fixed
  `log2(x + 1)` transform plus per-sample median centering). `"Log2Median"`
  is unchanged and still available; `"VSN"` now performs real VSN
  calibration. Documentation across README, the vignette, the Shiny method
  reference and the Rd pages has been corrected accordingly.
- `normalize_vsn()` validates its input before calibrating and never repairs
  data silently: infinite values, all-missing samples or features, fewer than
  2 samples, fewer than 42 lipid features, and non-numeric input all raise
  descriptive errors; missing values are reported and returned as `NA`
  (never imputed); negative values are passed to vsn unchanged but warn.
- If **vsn** is not installed, `"VSN"` raises an actionable error pointing to
  `BiocManager::install('vsn')`. It never falls back to another
  normalization method, and the Shiny app surfaces the error instead of
  reporting that normalization succeeded.
- The Shiny app now names the methods that were actually applied in its
  success notification and in the generated report, and the pipeline
  comparison commits nothing if either pipeline fails.

# LIPIDIFy 0.99.0

## New features
- Added `normalize_mean()` and full pipeline builder supporting TIC, PQN,
  Quantile, Log2Median, Median, Mean, Log2, Log10, Sqrt, None
- Added `get_normalization_descriptions()` for UI help text
- Shiny app: Select All / Deselect All buttons on every checkbox group
- Shiny app: group-coloured normalization comparison plots
- Shiny app: full plot history tracked across the session for complete reports
- Shiny app: samples sorted by group in Lipid Expression plots
- Shiny app: descriptive, timestamped download filenames for all plots

## Bug fixes
- `create_pca_plot_with_ellipses()`: fixed duplicate legend keys when ellipses
  are drawn (fill + color legends now merged into one)
- `create_lipid_expression_barplot()`: samples now sorted by group, not by
  sample name
- Excel export: worksheet names truncated to Excel's 31-character limit
- Normalization help text modal and ellipse help modal now display correctly
- `perform_plsda()`: replaced `cat()` with `message()` (Bioconductor requirement)
- `perform_enrichment_analysis()`: replaced `cat()` with `message()`

## Notes
- Quantile normalization producing near-identical boxplots is expected and
  correct behaviour; inline documentation and UI now make this explicit
- Median and Mean normalization can look similar on symmetric (log-normal)
  lipidomics data; on skewed data they differ — this is also expected
