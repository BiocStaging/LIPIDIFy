# Helper: truncate a string to n characters (for Excel sheet names, max 31)
.truncate_sheet_name <- function(x, n = 31) {
  x_clean <- gsub("[^A-Za-z0-9_.-]", "_", x)
  if (nchar(x_clean) > n) substr(x_clean, 1, n) else x_clean
}

# Helper: make a timestamped download filename
.dl_name <- function(..., ext = "png") {
  parts <- c(..., format(Sys.Date(), "%Y%m%d"))
  paste0(paste(parts, collapse = "_"), ".", ext)
}

# ---------------------------------------------------------------------------
# UI helper: checkbox group with Select All / Deselect All buttons
# ---------------------------------------------------------------------------
#' Checkbox Group with Select/Deselect All Buttons
#'
#' @param inputId  The checkbox group input ID.
#' @param label    Label above the checkbox group.
#' @param choices  Named character vector of choices (passed to
#'   \code{checkboxGroupInput}).
#' @param selected Initially selected choices.
#' @return A \code{tagList} containing the buttons and the checkbox group.
#' @noRd
.checkbox_group_with_buttons <- function(inputId, label, choices,
                                        selected = NULL) {
  btn_select_id <- paste0(inputId, "_select_all")
  btn_deselect_id <- paste0(inputId, "_deselect_all")

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::actionButton(btn_select_id,
        "Select All",
        class = "btn btn-xs btn-default",
        style = "margin-bottom:4px"
      )),
      shiny::column(6, shiny::actionButton(btn_deselect_id,
        "Deselect All",
        class = "btn btn-xs btn-default",
        style = "margin-bottom:4px"
      ))
    ),
    shiny::checkboxGroupInput(inputId, label,
      choices  = choices,
      selected = selected
    )
  )
}

# ---------------------------------------------------------------------------
# Server helper: wire Select All / Deselect All observers
# ---------------------------------------------------------------------------
#' Register Select-All and Deselect-All Observers for a CheckboxGroupInput
#'
#' @param session Shiny session object.
#' @param input   Shiny input object.
#' @param inputId The checkbox group input ID (buttons follow naming convention
#'   \code{<inputId>_select_all} and \code{<inputId>_deselect_all}).
#' @param choices_reactive A reactive expression (or function) that returns the
#'   current set of choices.
#' @noRd
.register_select_all_observers <- function(session, input, inputId,
                                          choices_reactive) {
  btn_select_id <- paste0(inputId, "_select_all")
  btn_deselect_id <- paste0(inputId, "_deselect_all")

  shiny::observeEvent(input[[btn_select_id]],
    {
      shiny::updateCheckboxGroupInput(session, inputId, selected = choices_reactive())
    },
    ignoreNULL = TRUE
  )

  shiny::observeEvent(input[[btn_deselect_id]],
    {
      shiny::updateCheckboxGroupInput(session, inputId, selected = character(0))
    },
    ignoreNULL = TRUE
  )
}

# ===========================================================================
# Main app function
# ===========================================================================

#' Launch the LIPIDIFy Shiny Application
#'
#' Starts an interactive Shiny dashboard for end-to-end lipidomics analysis,
#' including data upload, lipid classification, normalization, differential
#' analysis, enrichment analysis, and report generation.
#'
#' @param port Integer. Port number for the Shiny server. \code{NULL} lets
#'   Shiny pick a free port automatically.
#' @return Launches the Shiny application (does not return a value).
#' @export
#' @examples
#' if (interactive()) {
#'   launch_lipidomics_app()
#' }
launch_lipidomics_app <- function(port = NULL) {
  # ==========================================================================
  # UI
  # ==========================================================================
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "LIPIDIFy"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Welcome", tabName = "welcome", icon = shiny::icon("info-circle")),
        shinydashboard::menuItem("Data Upload", tabName = "upload", icon = shiny::icon("upload")),
        shinydashboard::menuItem("Lipid Classification", tabName = "classification", icon = shiny::icon("tags")),
        shinydashboard::menuItem("Raw Data Visualization", tabName = "raw_viz", icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Normalization", tabName = "normalization", icon = shiny::icon("balance-scale")),
        shinydashboard::menuItem("Preprocessing", tabName = "preprocessing", icon = shiny::icon("filter")),
        shinydashboard::menuItem("Normalized Data Visualization", tabName = "norm_viz", icon = shiny::icon("chart-bar")),
        shinydashboard::menuItem("Lipid Expression", tabName = "lipid_expression", icon = shiny::icon("flask")),
        shinydashboard::menuItem("Differential Analysis", tabName = "diff_analysis", icon = shiny::icon("calculator")),
        shinydashboard::menuItem("Results Visualization", tabName = "results_viz", icon = shiny::icon("chart-area")),
        shinydashboard::menuItem("Enrichment Analysis", tabName = "enrichment", icon = shiny::icon("search-plus")),
        shinydashboard::menuItem("Enrichment Visualization", tabName = "enrichment_viz", icon = shiny::icon("network-wired")),
        shinydashboard::menuItem("Generate Report", tabName = "report", icon = shiny::icon("file-alt"))
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        # ------------------------------------------------------------------
        # Welcome
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "welcome",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Welcome to LIPIDIFy", status = "primary",
              solidHeader = TRUE, width = 12,
              shiny::HTML('
<style>
/* ---- welcome page styles ---- */
.lf-section  { margin-bottom: 28px; }
.lf-section h3 { color: #2c3e50; margin-top: 0; }
.lf-section h4 { color: #2980b9; border-bottom: 2px solid #2980b9;
                 padding-bottom: 5px; margin-top: 20px; }
.lf-section h5 { color: #27ae60; margin: 12px 0 4px 0; font-size: 14px; }
.lf-panel    { background: #f8f9fa; border-left: 4px solid #2980b9;
               padding: 12px 16px; border-radius: 4px; margin: 8px 0 14px 0; }
.lf-panel.green  { border-color: #27ae60; }
.lf-panel.orange { border-color: #e67e22; }
.lf-panel.red    { border-color: #e74c3c; }
.lf-code     { font-family: "Courier New", monospace; font-size: 11.5px;
               background: #2c3e50; color: #ecf0f1; padding: 12px 14px;
               border-radius: 4px; overflow-x: auto; margin: 8px 0;
               white-space: pre; }
.lf-tip      { background: #d5f5e3; border: 1px solid #a9dfbf;
               border-radius: 4px; padding: 10px 14px; margin: 10px 0; }
.lf-warn     { background: #fef9e7; border: 1px solid #f9e79f;
               border-radius: 4px; padding: 10px 14px; margin: 10px 0; }
.lf-danger   { background: #fdedec; border: 1px solid #f1948a;
               border-radius: 4px; padding: 10px 14px; margin: 10px 0; }
.lf-badge    { display: inline-block; background: #2980b9; color: white;
               border-radius: 12px; padding: 1px 9px; font-size: 11px;
               font-weight: bold; margin-right: 4px; }
.lf-badge.green  { background: #27ae60; }
.lf-badge.orange { background: #e67e22; }
.lf-badge.red    { background: #e74c3c; }
table.lf-table { width: 100%; border-collapse: collapse; font-size: 12.5px; margin: 8px 0; }
table.lf-table th { background: #2980b9; color: white; padding: 7px 10px; text-align: left; }
table.lf-table td { padding: 6px 10px; border-bottom: 1px solid #ddd; vertical-align: top; }
table.lf-table tr:nth-child(even) td { background: #f2f3f4; }
</style>

<!-- ===================== HEADER ===================== -->
<div class="lf-section">
<h3><i class="fa fa-flask"></i> LIPIDIFy - Complete Lipidomics Analysis Platform</h3>
<p style="font-size:14px">
LIPIDIFy is an interactive R/Shiny application covering the full lipidomics workflow:
from raw data upload through normalisation, differential abundance testing, enrichment
analysis, and publication-ready report generation.
No programming knowledge is required.
</p>
<div class="lf-tip">
<i class="fa fa-play-circle"></i>
<strong>New user? Start here:</strong> Click <em>Data Upload</em> in the left menu,
then click <strong>"Load Example Data"</strong>. This loads a built-in demo dataset
so you can explore every feature of the app before uploading your own files.
</div>
</div>

<!-- ===================== WORKFLOW ===================== -->
<div class="lf-section">
<h4><i class="fa fa-list-ol"></i> Step-by-Step Workflow</h4>

<table class="lf-table">
<tr><th style="width:5%">#</th><th style="width:22%">Section</th><th>What it does</th><th style="width:28%">What you need</th></tr>
<tr>
  <td><span class="lf-badge">1</span></td>
  <td><strong>Data Upload</strong></td>
  <td>Load your lipidomics CSV. The app separates metadata (sample names, groups) from
      lipid abundance columns automatically.</td>
  <td>A CSV file - see format below</td>
</tr>
<tr>
  <td><span class="lf-badge">2</span></td>
  <td><strong>Lipid Classification</strong></td>
  <td>Each lipid is classified by class, subclass, and fatty-acid saturation (SFA/MUFA/PUFA)
      from its name. You can upload your own classification table or reset to automatic
      at any time.</td>
  <td>Optional: custom classification CSV</td>
</tr>
<tr>
  <td><span class="lf-badge">3</span></td>
  <td><strong>Raw Data Visualization</strong></td>
  <td>Inspect intensity distributions before any processing. Filter by group,
      choose boxplot / violin / density / histogram. Helps spot batch effects
      and outlier samples.</td>
  <td>Data loaded (step 1)</td>
</tr>
<tr>
  <td><span class="lf-badge">4</span></td>
  <td><strong>Normalisation</strong></td>
  <td>Build two normalisation pipelines (chain methods in any order) and compare them
      side-by-side with boxplot, violin, or density plots. Apply the better pipeline
      before downstream analysis.</td>
  <td>Data loaded (step 1)</td>
</tr>
<tr>
  <td><span class="lf-badge">5</span></td>
  <td><strong>Normalised Data Visualization</strong></td>
  <td>Verify normalisation quality with boxplots, PCA, or PLS-DA. Filter groups,
      add confidence ellipses, show sample name labels on PCA.</td>
  <td>Normalisation applied (step 4)</td>
</tr>
<tr>
  <td><span class="lf-badge">6</span></td>
  <td><strong>Lipid Expression</strong></td>
  <td>Bar plots of individual lipid abundance across samples, automatically sorted
      by group. Select up to 10 lipids at once.</td>
  <td>Data loaded (step 1)</td>
</tr>
<tr>
  <td><span class="lf-badge">7</span></td>
  <td><strong>Differential Analysis</strong></td>
  <td>Find significantly changed lipids between groups using limma or EdgeR.
      All pairwise contrasts are created automatically. Results include
      log-fold change and FDR-corrected p-values.</td>
  <td>Normalisation applied (step 4)</td>
</tr>
<tr>
  <td><span class="lf-badge">8</span></td>
  <td><strong>Results Visualization</strong></td>
  <td>Volcano plots coloured by lipid class, or heatmaps of the most significant
      features. Download as PNG or PDF.</td>
  <td>Differential analysis run (step 7)</td>
</tr>
<tr>
  <td><span class="lf-badge">9</span></td>
  <td><strong>Enrichment Analysis</strong></td>
  <td>GSEA-style enrichment on lipid classes using fgsea. Visualise as dot plots
      or bar plots. Upload your own biological sets for custom enrichment.</td>
  <td>Differential analysis run (step 7)</td>
</tr>
<tr>
  <td><span class="lf-badge">10</span></td>
  <td><strong>Generate Report</strong></td>
  <td>Export a complete HTML or PDF report with all figures, summary tables,
      and session information. Every plot generated during the session is included.</td>
  <td>Any analyses completed</td>
</tr>
</table>
</div>

<!-- ===================== FILE 1: MAIN DATA ===================== -->
<div class="lf-section">
<h4><i class="fa fa-file-csv"></i> Input File 1 - Main Lipidomics Data (required)</h4>

<p>Format: <strong>CSV</strong> (comma-separated values). Save from Excel as
"CSV UTF-8" or "CSV (Comma delimited)".</p>

<div class="lf-panel">
<strong>Rules:</strong>
<ul style="margin: 6px 0 0 0">
  <li>Each <strong>row</strong> = one sample</li>
  <li>First columns = <strong>metadata</strong> (Sample Name, Sample Group, etc.)</li>
  <li>Remaining columns = <strong>lipid abundance values</strong> (numbers)</li>
  <li>"Sample Name" and "Sample Group" are the two most important metadata columns</li>
  <li>All lipid abundance values must be positive numbers (no text, no % signs)</li>
  <li>Missing values should be left empty or written as NA</li>
</ul>
</div>

<h5>Minimal example (2 metadata columns + 3 lipids):</h5>
<div class="lf-code">Sample Name,Sample Group,PC 16:0_18:1,PE 18:0_20:4,TG 16:0_18:1_20:4
Control_1,Control,15420.5,8932.1,45231.0
Control_2,Control,14856.2,9125.4,43870.3
Control_3,Control,16201.8,8745.2,44502.7
Treated_1,Treatment,18932.1,7845.6,38921.4
Treated_2,Treatment,19201.4,7612.3,37845.2
Treated_3,Treatment,17890.6,8012.5,39234.1</div>

<h5>Extended example (4 metadata columns):</h5>
<div class="lf-code">Sample Name,Sample Group,Tumour ID,Weight (mg),PC 16:0_18:1,PE 18:0_20:4,SM 18:1
P001_pre,Pre-treatment,T001,52.3,15420.5,8932.1,12450.2
P001_post,Post-treatment,T001,48.7,19832.4,7541.2,10231.5
P002_pre,Pre-treatment,T002,55.1,14756.3,9012.4,13102.8
P002_post,Post-treatment,T002,51.2,20145.7,7234.8,9845.3</div>

<div class="lf-tip">
<i class="fa fa-lightbulb-o"></i>
<strong>Tip - Lipid naming:</strong> Use standard lipid notation so the automatic
classification works correctly. Examples:
</div>
<table class="lf-table" style="margin-top:4px">
<tr><th>Class</th><th>Example name</th><th>Meaning</th></tr>
<tr><td>Phosphatidylcholine</td><td>PC 16:0_18:1</td><td>PC with 16:0 and 18:1 acyl chains</td></tr>
<tr><td>Phosphatidylethanolamine</td><td>PE 18:0_20:4</td><td>PE with 18:0 and 20:4 chains</td></tr>
<tr><td>Triacylglycerol</td><td>TG 16:0_18:1_20:4</td><td>TG with three chains</td></tr>
<tr><td>Lysophosphatidylcholine</td><td>LPC 18:2</td><td>LPC with one 18:2 chain</td></tr>
<tr><td>Sphingomyelin</td><td>SM 18:1</td><td>SM with 18:1 chain</td></tr>
<tr><td>Ceramide</td><td>Cer 16:0</td><td>Ceramide with 16:0 chain</td></tr>
<tr><td>Cholesteryl ester</td><td>CE 18:1</td><td>CE with 18:1 chain</td></tr>
<tr><td>Diacylglycerol</td><td>DG 16:0_18:1</td><td>DG with two chains</td></tr>
</table>

<div class="lf-warn" style="margin-top:10px">
<i class="fa fa-exclamation-triangle"></i>
<strong>Common mistakes:</strong>
<ul style="margin: 4px 0 0 0">
  <li>Do not include a "Total" or "Sum" row at the bottom</li>
  <li>Do not include QC pool rows (rows with "PBQC" in Sample Group are removed automatically)</li>
  <li>Do not mix text and numbers in the same lipid column</li>
  <li>Sample names must be unique - no two rows with the same Sample Name</li>
</ul>
</div>
</div>

<!-- ===================== FILE 2: CLASSIFICATION ===================== -->
<div class="lf-section">
<h4><i class="fa fa-tags"></i> Input File 2 - Custom Lipid Classification (optional)</h4>

<p>By default, the app classifies lipids automatically from their names. If you have
your own classification (e.g., from LipidMaps, or a custom biological grouping), you
can upload a CSV to override it.</p>

<div class="lf-panel green">
<strong>Rules:</strong>
<ul style="margin: 6px 0 0 0">
  <li>Format: <strong>CSV</strong></li>
  <li>First column must be named exactly <code>Lipid</code></li>
  <li>Lipid names must match your data file <strong>exactly</strong> (case-sensitive, character-for-character)</li>
  <li>You can have <strong>any number of additional columns</strong> - each becomes a coloring option in volcano plots and an enrichment category</li>
  <li>Each lipid should appear only once (duplicates are removed automatically, keeping the first row)</li>
</ul>
</div>

<h5>Simple example (one classification column):</h5>
<div class="lf-code">Lipid,LipidClass
LPC(13:0),LPC Lipid Class
LPC(14:0),LPC Lipid Class
LPC(16:0),LPC Lipid Class
PC 16:0_18:1,Phospholipid
PE 18:0_20:4,Phospholipid
TG 16:0_18:1_20:4,Glycerolipid
Cer 16:0,Sphingolipid</div>

<h5>Extended example (multiple classification columns):</h5>
<div class="lf-code">Lipid,LipidClass,BiologicalRole,Membrane
PC 16:0_18:1,Phospholipid,Structural,Yes
PE 18:0_20:4,Phospholipid,Signaling,Yes
TG 16:0_18:1_20:4,Glycerolipid,Energy Storage,No
Cer 16:0,Sphingolipid,Apoptosis,No
SM 18:1,Sphingolipid,Structural,Yes</div>

<div class="lf-danger">
<i class="fa fa-exclamation-circle"></i>
<strong>Critical - exact name matching:</strong> If your classification says <code>LPC(16:0)</code>
but your data column is named <code>LPC 16:0</code>, no match will occur and all lipids will
appear as "Unknown" in the volcano plot. The app will show you a diagnostic message with
example names from both sides to help you find the mismatch.
</div>

<div class="lf-tip">
<i class="fa fa-undo"></i>
<strong>Reset:</strong> After loading a custom classification, you can always click
"Reset to Automatic Classification" to go back to the default name-based classification.
</div>
</div>

<!-- ===================== FILE 3: ENRICHMENT SETS ===================== -->
<div class="lf-section">
<h4><i class="fa fa-search-plus"></i> Input File 3 - Custom Enrichment Sets (optional)</h4>

<p>Enrichment analysis tests whether lipids belonging to a particular biological set
(pathway, class, function) are collectively increased or decreased. By default, the app
uses the lipid classification columns for enrichment. You can also upload your own
biologically meaningful sets.</p>

<div class="lf-panel orange">
<strong>Rules:</strong>
<ul style="margin: 6px 0 0 0">
  <li>Format: <strong>CSV</strong></li>
  <li>Exactly two columns: <code>Lipid</code> and <code>Set_Name</code></li>
  <li>One row per lipid-per-set assignment</li>
  <li>A lipid can appear in <strong>multiple rows</strong> to belong to multiple sets - this is normal and expected</li>
  <li>Set names can be anything descriptive (pathway names, functional categories, etc.)</li>
  <li>Lipid names must match your data file exactly (same rule as classification)</li>
  <li>Minimum set size is 3 lipids by default (configurable in the Enrichment tab)</li>
</ul>
</div>

<h5>Example - pathway-based sets:</h5>
<div class="lf-code">Lipid,Set_Name
PC 16:0_18:1,Membrane_Phospholipids
PC 18:0_20:4,Membrane_Phospholipids
PE 18:0_20:4,Membrane_Phospholipids
PE 16:0_18:1,Membrane_Phospholipids
TG 16:0_18:1_18:1,Energy_Storage
TG 16:0_18:0_18:1,Energy_Storage
TG 18:0_18:1_18:1,Energy_Storage
Cer 16:0,Apoptosis_Ceramides
Cer 18:0,Apoptosis_Ceramides
Cer 24:1,Apoptosis_Ceramides
PE 18:0_20:4,AA_Derived_Lipids
PC 18:0_20:4,AA_Derived_Lipids</div>

<p style="font-size:12px; color:#666; margin-top:4px">
Note: PE 18:0_20:4 appears in both "Membrane_Phospholipids" and "AA_Derived_Lipids" - this
is correct. A lipid can belong to multiple biological sets simultaneously.
</p>

<h5>Example - saturation-based sets:</h5>
<div class="lf-code">Lipid,Set_Name
PC 18:0_18:0,Saturated_Species
PE 18:0_18:0,Saturated_Species
SM 16:0,Saturated_Species
PC 16:0_18:1,Monounsaturated_Species
PE 18:0_18:1,Monounsaturated_Species
PC 18:0_20:4,Polyunsaturated_Species
PE 18:0_20:4,Polyunsaturated_Species
PC 18:0_22:6,Polyunsaturated_Species</div>
</div>

<!-- ===================== NORMALISATION REFERENCE ===================== -->
<div class="lf-section">
<h4><i class="fa fa-cogs"></i> Normalisation Methods Reference</h4>

<table class="lf-table">
<tr><th style="width:15%">Method</th><th style="width:42%">What it does</th><th>When to use it</th><th style="width:12%">Common use</th></tr>
<tr>
  <td><strong>TIC</strong></td>
  <td>Divides each sample by its total signal (sum of all lipids), then rescales to the global mean. Compensates for differences in total injected material.</td>
  <td>When injection volumes or sample amounts vary between samples</td>
  <td><span class="lf-badge">Common</span></td>
</tr>
<tr>
  <td><strong>PQN</strong></td>
  <td>Probabilistic Quotient Normalization. Uses the median spectrum as reference and normalises by the median quotient. More robust than TIC when some lipids change dramatically.</td>
  <td>General-purpose; preferred when large fold-changes are expected</td>
  <td><span class="lf-badge">Common</span></td>
</tr>
<tr>
  <td><strong>Quantile</strong></td>
  <td>Forces all samples to have identical intensity distributions. After this, all boxplots will look the same - <em>this is the expected, correct behaviour</em>, not a bug.</td>
  <td>When you want to remove all between-sample distribution differences</td>
  <td><span class="lf-badge orange">Advanced</span></td>
</tr>
<tr>
  <td><strong>VSN</strong></td>
  <td>Variance Stabilizing Normalization (simplified): log2 transform + median centering. Stabilises variance that scales with intensity.</td>
  <td>When low-abundance lipids have high variance</td>
  <td><span class="lf-badge orange">Advanced</span></td>
</tr>
<tr>
  <td><strong>Median</strong></td>
  <td>Scales each sample so its median equals the global median across all samples. Robust to outlier lipids.</td>
  <td>Simple general-purpose normalisation, good first choice</td>
  <td><span class="lf-badge">Common</span></td>
</tr>
<tr>
  <td><strong>Mean</strong></td>
  <td>Same as Median but uses the mean. Slightly more sensitive to outliers. On symmetric data, results look similar to Median - this is expected.</td>
  <td>When you prefer mean-based scaling</td>
  <td><span class="lf-badge">Common</span></td>
</tr>
<tr>
  <td><strong>Log2</strong></td>
  <td>Log base-2 transformation. Reduces the right skew typical of MS intensity data. Makes the data more symmetrically distributed.</td>
  <td>Almost always applied as the last step in a pipeline</td>
  <td><span class="lf-badge green">Recommended</span></td>
</tr>
<tr>
  <td><strong>Log10</strong></td>
  <td>Log base-10 transformation. Same idea as Log2 but gives a different scale. Useful when reporting fold-changes in log10 units.</td>
  <td>Alternative to Log2</td>
  <td><span class="lf-badge">Optional</span></td>
</tr>
<tr>
  <td><strong>Sqrt</strong></td>
  <td>Square-root transformation. A milder alternative to log for data that is less skewed.</td>
  <td>Count-like data or when log produces too-compressed values</td>
  <td><span class="lf-badge">Optional</span></td>
</tr>
<tr>
  <td><strong>None</strong></td>
  <td>No transformation applied. Passes data through unchanged.</td>
  <td>When data is already normalised externally</td>
  <td></td>
</tr>
</table>

<div class="lf-panel" style="margin-top:10px">
<strong>Recommended pipelines for beginners:</strong>
<ul style="margin: 6px 0 0 0">
  <li><strong>TIC + Log2</strong> - most common starting point for MS lipidomics</li>
  <li><strong>PQN + Log2</strong> - robust alternative, especially with heterogeneous samples</li>
  <li><strong>Median + Log2</strong> - simple and effective when injection amounts are consistent</li>
</ul>
Methods are applied left-to-right in the order you select them.
</div>
</div>

<!-- ===================== STATISTICAL ANALYSIS REFERENCE ===================== -->
<div class="lf-section">
<h4><i class="fa fa-calculator"></i> Statistical Analysis Reference</h4>

<table class="lf-table">
<tr><th style="width:18%">Term</th><th>Definition</th></tr>
<tr>
  <td><strong>logFC</strong></td>
  <td>Log2 Fold Change. Positive = higher in the first group; negative = lower.
      A logFC of 1 means 2-fold higher; logFC of 2 means 4-fold higher; logFC of -1 means 2-fold lower.</td>
</tr>
<tr>
  <td><strong>adj.P.Val</strong></td>
  <td>Adjusted p-value (False Discovery Rate, Benjamini-Hochberg method).
      Corrects for testing many lipids simultaneously. Use this, not the raw p-value.
      Threshold: typically 0.05 (5% FDR).</td>
</tr>
<tr>
  <td><strong>Contrast</strong></td>
  <td>A pairwise comparison, written as "GroupA - GroupB". Positive logFC means
      higher in GroupA. Example: "Treatment - Control" means the treatment group
      is compared to control.</td>
</tr>
<tr>
  <td><strong>NES</strong></td>
  <td>Normalized Enrichment Score (enrichment analysis). Positive = set enriched in
      increased lipids; negative = set enriched in decreased lipids.</td>
</tr>
<tr>
  <td><strong>padj (enrichment)</strong></td>
  <td>Adjusted p-value for enrichment results. A pathway is considered significantly
      enriched when padj &lt; 0.25 (standard GSEA threshold).</td>
</tr>
</table>
</div>
              ')
            )
          )
        ),
        # ------------------------------------------------------------------
        # Data Upload
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "upload",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Upload Lipidomics Data", status = "primary",
              solidHeader = TRUE, width = 6,
              shiny::fileInput("file", "Choose CSV File", accept = ".csv"),
              .checkbox_group_with_buttons(
                "metadata_cols", "Metadata Columns:",
                choices = c("Sample Name", "Sample Group", "Tumour ID", "Weight (mg)"),
                selected = c("Sample Name", "Sample Group")
              ),
              shiny::actionButton("load_data", "Load Data", class = "btn-primary"),
              shiny::br(), shiny::br(),
              shiny::actionButton("load_example", "Load Example Data", class = "btn-info")
            ),
            shinydashboard::box(
              title = "Data Summary", status = "info",
              solidHeader = TRUE, width = 6,
              shiny::verbatimTextOutput("data_summary"),
              DT::dataTableOutput("data_preview")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Classification
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "classification",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Lipid Classification", status = "primary",
              solidHeader = TRUE, width = 6,
              shiny::h4("Automatic Classification"),
              shiny::p("Lipids are classified automatically from their names."),
              shiny::actionLink("info_classification", " Help",
                icon = shiny::icon("info-circle")
              ),
              shiny::br(), shiny::br(),
              shiny::downloadButton("download_classification",
                "Download Classification",
                class = "btn-success"
              ),
              shiny::hr(),
              shiny::h4("Upload Custom Classification"),
              shiny::fileInput("custom_classification_file",
                "Choose CSV File",
                accept = ".csv"
              ),
              shiny::helpText("File must have a 'Lipid' column plus classification columns."),
              shiny::actionButton("load_custom_classification",
                "Load Custom Classification",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::actionButton("reset_classification",
                "Reset to Automatic Classification",
                class = "btn-warning",
                icon  = shiny::icon("undo")
              )
            ),
            shinydashboard::box(
              title = "Current Classification", status = "info",
              solidHeader = TRUE, width = 6,
              DT::dataTableOutput("classification_table")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Raw Data Visualization
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "raw_viz",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Raw Data Visualization Options", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::radioButtons("view_mode", "View Mode:",
                choices = c(
                  "By Sample" = "sample",
                  "By Lipid" = "lipid"
                ),
                selected = "sample"
              ),
              shiny::selectInput("plot_type", "Plot Type:",
                choices = c(
                  "Boxplot"   = "boxplot",
                  "Violin"    = "violin",
                  "Density"   = "density",
                  "Histogram" = "histogram",
                  "Heatmap"   = "heatmap"
                )
              ),
              shiny::conditionalPanel(
                condition = "input.plot_type == 'heatmap'",
                shiny::numericInput("raw_heatmap_top_n",
                  "Top N Variable Lipids:",
                  value = 50, min = 10, max = 200
                )
              ),
              shiny::conditionalPanel(
                condition = "input.view_mode == 'lipid'",
                shiny::numericInput("top_n_lipids", "Top N Variable Lipids:",
                  value = 30, min = 5, max = 200
                )
              ),
              shiny::conditionalPanel(
                condition = "input.view_mode == 'sample'",
                shiny::checkboxInput("limit_top_samples",
                  "Limit to top N samples", value = FALSE
                ),
                shiny::conditionalPanel(
                  condition = "input.limit_top_samples == true",
                  shiny::numericInput("top_n_samples", "Top N Variable Samples:",
                    value = 30, min = 2, max = 1000
                  )
                )
              ),
              shiny::hr(),
              shiny::h5("Filter and Colour:"),
              .checkbox_group_with_buttons(
                "raw_filter_groups", "Show Groups Only:",
                choices = NULL, selected = NULL
              ),
              shiny::checkboxInput("raw_color_by_group",
                "Colour by Sample Group",
                value = TRUE
              ),
              shiny::actionButton("create_raw_plot", "Create Plot",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::selectInput(
                "img_format_raw", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_raw_plot", "Download Plot",
                class = "btn-success"
              )
            ),
            shinydashboard::box(
              title = "Raw Data Visualization", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::uiOutput("raw_plot_ui")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Normalization
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "normalization",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Normalization Pipeline Builder", status = "primary",
              solidHeader = TRUE, width = 4,

              # Help button
              shiny::actionLink("info_norm_methods", " Normalization Method Descriptions",
                icon = shiny::icon("info-circle")
              ),
              shiny::hr(),
              shiny::h4("Pipeline 1:"),
              .checkbox_group_with_buttons(
                "norm_methods_1", "Select Methods (applied in order):",
                choices = get_normalization_methods(),
                selected = c("TIC", "Log2")
              ),
              shiny::h4("Pipeline 2 (for comparison):"),
              .checkbox_group_with_buttons(
                "norm_methods_2", "Select Methods (applied in order):",
                choices = get_normalization_methods(),
                selected = "PQN"
              ),
              shiny::hr(),
              shiny::h4("Comparison Plot Options"),
              shiny::selectInput("norm_compare_plot_type",
                "Plot Type:",
                choices = c(
                  "Boxplot" = "boxplot",
                  "Violin" = "violin",
                  "Density" = "density"
                )
              ),
              shiny::conditionalPanel(
                condition = "input.norm_compare_plot_type != 'density'",
                shiny::radioButtons("norm_compare_view_mode", "View By:",
                  choices = c("By Sample" = "sample", "By Lipid" = "lipid"),
                  selected = "sample", inline = TRUE
                ),
                shiny::checkboxInput("norm_compare_limit_top",
                  "Limit to top N (most variable)", value = FALSE
                ),
                shiny::conditionalPanel(
                  condition = "input.norm_compare_limit_top == true",
                  shiny::numericInput("norm_compare_top_n", "Top N:",
                    value = 30, min = 2, max = 1000
                  )
                )
              ),
              shiny::checkboxInput("norm_color_by_group",
                "Colour boxes by Sample Group",
                value = TRUE
              ),
              shiny::actionButton("compare_pipelines", "Compare Pipelines",
                class = "btn-info"
              ),
              shiny::hr(),
              shiny::selectInput("chosen_pipeline",
                "Apply Pipeline:",
                choices = c(
                  "Pipeline 1" = "1",
                  "Pipeline 2" = "2"
                )
              ),
              shiny::actionButton("apply_normalization",
                "Apply Selected Pipeline",
                class = "btn-primary"
              )
            ),
            shinydashboard::box(
              title = "Pipeline Comparison", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::plotOutput("pipeline_comparison", height = "650px"),
              shiny::br(),
              shiny::selectInput(
                "img_format_pipeline", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_pipeline_comparison",
                "Download Comparison",
                class = "btn-success"
              )
            )
          )
        ),

        # ------------------------------------------------------------------
        # Preprocessing (Imputation + Batch correction) - applied post-normalisation
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "preprocessing",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Missing Value Imputation", status = "primary",
              solidHeader = TRUE, width = 6,
              shiny::p(
                "Replace missing values (NA) in the ", shiny::strong("normalised"),
                " data. Imputing after normalisation ensures imputed values ",
                "are on the same scale as observed values."
              ),
              shiny::verbatimTextOutput("missing_value_summary"),
              shiny::hr(),
              shiny::selectInput(
                "imputation_method", "Imputation Method:",
                choices = get_imputation_methods(),
                selected = "half_min"
              ),
              shiny::conditionalPanel(
                condition = "input.imputation_method == 'knn'",
                shiny::numericInput("imputation_k",
                  "Number of Neighbours (k):",
                  value = 5L, min = 2L, max = 20L
                )
              ),
              shiny::actionLink("info_imputation", " Method Descriptions",
                icon = shiny::icon("info-circle")
              ),
              shiny::br(), shiny::br(),
              shiny::actionButton("run_imputation", "Apply Imputation",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::actionButton("reset_imputation",
                "Reset to Pre-Imputation Data",
                class  = "btn-warning",
                icon   = shiny::icon("undo")
              )
            ),
            shinydashboard::box(
              title = "Batch Effect Correction", status = "primary",
              solidHeader = TRUE, width = 6,
              shiny::p(
                "Remove known technical batch effects while preserving ",
                "biological signal."
              ),
              shiny::helpText(
                "Requires a dedicated batch column in your metadata ",
                "(e.g. 'Run', 'Plate', 'Batch'). If your CSV only has ",
                "'Sample Name' and 'Sample Group', add a batch column ",
                "and reload the data before using this feature."
              ),
              shiny::selectInput("batch_column",
                "Batch Column in Metadata:", choices = NULL
              ),
              shiny::selectInput("batch_method", "Method:",
                choices = c(
                  "limma - removeBatchEffect (recommended)" = "limma",
                  "ComBat - sva package (robust for large effects)"  = "combat"
                )
              ),
              shiny::selectInput("batch_group_column",
                "Group Column to Protect:", choices = NULL
              ),
              shiny::actionLink("info_batch", " Batch Correction Help",
                icon = shiny::icon("info-circle")
              ),
              shiny::br(), shiny::br(),
              shiny::actionButton("run_batch_correction",
                "Apply Batch Correction",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::actionButton("reset_batch_correction",
                "Reset to Pre-Batch Data",
                class = "btn-warning",
                icon  = shiny::icon("undo")
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Preprocessing Status", status = "info",
              solidHeader = TRUE, width = 12,
              shiny::verbatimTextOutput("preprocessing_status")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Normalized Data Visualization
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "norm_viz",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Normalized Data Visualization", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectInput("norm_plot_type", "Plot Type:",
                choices = c(
                  "Boxplot" = "boxplot",
                  "Violin"  = "violin",
                  "Density" = "density",
                  "Heatmap" = "heatmap",
                  "PCA"     = "pca",
                  "PLS-DA"  = "plsda"
                )
              ),
              shiny::conditionalPanel(
                condition = "input.norm_plot_type == 'heatmap'",
                shiny::numericInput("norm_heatmap_top_n",
                  "Top N Variable Lipids:",
                  value = 50, min = 10, max = 200
                )
              ),
              shiny::conditionalPanel(
                condition = paste0(
                  "input.norm_plot_type == 'boxplot' || ",
                  "input.norm_plot_type == 'violin' || ",
                  "input.norm_plot_type == 'density'"
                ),
                shiny::checkboxInput("norm_limit_top_samples",
                  "Limit to top N samples", value = FALSE
                ),
                shiny::conditionalPanel(
                  condition = "input.norm_limit_top_samples == true",
                  shiny::numericInput("norm_top_n_samples",
                    "Top N Variable Samples:",
                    value = 30, min = 2, max = 1000
                  )
                )
              ),
              shiny::hr(),
              shiny::h5("Filter Samples:"),
              .checkbox_group_with_buttons(
                "norm_filter_groups", "Show Groups:",
                choices = NULL, selected = NULL
              ),
              shiny::checkboxInput("norm_color_by_group_viz",
                "Colour by Sample Group",
                value = TRUE
              ),
              shiny::conditionalPanel(
                condition = "input.norm_plot_type == 'pca' || input.norm_plot_type == 'plsda'",
                shiny::selectInput("group_column", "Group Column:", choices = NULL),
                shiny::actionLink("info_ellipses", " Ellipse Help",
                  icon = shiny::icon("info-circle")
                ),
                .checkbox_group_with_buttons(
                  "groups_included", "Include Groups:",
                  choices = NULL, selected = NULL
                ),
                shiny::radioButtons("ellipse_type", "Ellipse Type:",
                  choices = c(
                    "None" = "none",
                    "Confidence Ellipse" = "confidence",
                    "Visual Circle" = "visual"
                  ),
                  selected = "none"
                ),
                shiny::checkboxInput("show_sample_labels",
                  "Show Sample Labels",
                  value = FALSE
                )
              ),
              shiny::actionButton("create_norm_plot", "Create Plot",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::selectInput(
                "img_format_norm", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_norm_plot", "Download Plot",
                class = "btn-success"
              )
            ),
            shinydashboard::box(
              title = "Normalized Data Plot", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::uiOutput("norm_plot_ui")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Lipid Expression
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "lipid_expression",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Lipid Expression Settings", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectizeInput("selected_lipids", "Select Lipids:",
                choices = NULL, multiple = TRUE,
                options = list(
                  maxItems = 10,
                  placeholder = "Search for lipids..."
                )
              ),
              shiny::radioButtons("expression_selection_mode", "Select by:",
                choices = c(
                  "Samples" = "samples",
                  "Groups" = "groups"
                ),
                selected = "groups"
              ),
              shiny::conditionalPanel(
                condition = "input.expression_selection_mode == 'samples'",
                .checkbox_group_with_buttons(
                  "selected_samples", "Select Samples:",
                  choices = NULL, selected = NULL
                )
              ),
              shiny::conditionalPanel(
                condition = "input.expression_selection_mode == 'groups'",
                .checkbox_group_with_buttons(
                  "selected_expression_groups", "Select Groups:",
                  choices = NULL, selected = NULL
                )
              ),
              shiny::radioButtons("expression_data_type", "Data Type:",
                choices  = c("Raw" = "raw", "Normalized" = "normalized"),
                selected = "normalized"
              ),
              shiny::actionButton("create_expression_plot", "Create Plots",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::selectInput(
                "img_format_expression", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_expression_plots",
                "Download Plots",
                class = "btn-success"
              )
            ),
            shinydashboard::box(
              title = "Lipid Expression Plots", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::uiOutput("expression_plots_ui")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Differential Analysis
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "diff_analysis",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Differential Analysis Settings", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectInput("group_col_diff", "Group Column:", choices = NULL),
              shiny::hr(),
              shiny::h4("Analysis Method"),
              shiny::radioButtons("diff_method", "Select Method:",
                choices  = c("limma" = "limma", "EdgeR" = "edger"),
                selected = "limma"
              ),
              shiny::actionLink("info_methods", " Method Help",
                icon = shiny::icon("info-circle")
              ),
              shiny::hr(),
              shiny::h4("Select Contrasts"),
              shiny::actionLink("info_contrasts", " Contrast Help",
                icon = shiny::icon("info-circle")
              ),
              shiny::uiOutput("contrast_selection_ui"),
              shiny::hr(),
              shiny::h4("Custom Contrasts"),
              shiny::helpText(
                "Enter one contrast per line using limma notation. ",
                "Use this for complex comparisons such as collapsing groups:"
              ),
              shiny::tags$pre(
                style = "font-size:11px; background:#f5f5f5; padding:6px; border-radius:4px;",
                "(GroupA + GroupB)/2 - (GroupC + GroupD)/2"
              ),
              shiny::textAreaInput(
                "custom_contrasts_text",
                label    = NULL,
                value    = "",
                rows     = 3,
                placeholder = "e.g.\n(ND_Vehicle + ND_DabTram)/2 - (HFD_Vehicle + HFD_DabTram)/2"
              ),
              shiny::hr(),
              shiny::actionButton("run_diff_analysis", "Run Analysis",
                class = "btn-primary"
              )
            ),
            shinydashboard::box(
              title = "Analysis Summary", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::verbatimTextOutput("diff_summary"),
              shiny::hr(),
              shiny::selectInput("contrast_display",
                "Select Contrast to Display:",
                choices = NULL
              ),
              DT::dataTableOutput("diff_results_table"),
              shiny::br(),
              shiny::downloadButton("download_current_table",
                "Download Current Table (CSV)",
                class = "btn-info"
              ),
              shiny::downloadButton("download_results",
                "Download All Results (Excel)",
                class = "btn-success"
              )
            )
          )
        ),

        # ------------------------------------------------------------------
        # Results Visualization
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "results_viz",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Visualization Options", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectInput("contrast_select", "Select Contrast:", choices = NULL),
              shiny::selectInput("viz_type", "Visualization Type:",
                choices = c(
                  "Volcano Plot" = "volcano",
                  "Heatmap" = "heatmap"
                )
              ),
              shiny::conditionalPanel(
                condition = "input.viz_type == 'volcano'",
                shiny::numericInput("logfc_threshold", "LogFC Threshold:",
                  value = 1, min = 0, step = 0.1
                ),
                shiny::numericInput("pval_threshold", "P-value Threshold:",
                  value = 0.05, min = 0, max = 1, step = 0.01
                ),
                shiny::numericInput("top_labels", "Top N Labels:",
                  value = 15, min = 0, max = 50
                ),
                shiny::checkboxInput("color_by_class",
                  "Color Significant by Classification",
                  value = FALSE
                ),
                shiny::conditionalPanel(
                  condition = "input.color_by_class",
                  shiny::selectInput("color_column", "Color by:",
                    choices = c("LipidGroup", "LipidType", "Saturation")
                  )
                )
                # Note: choices are updated dynamically in the server when
                # classification data changes (handles custom classification columns)
              ),
              shiny::conditionalPanel(
                condition = "input.viz_type == 'heatmap'",
                shiny::numericInput("heatmap_top_n", "Top N Features:",
                  value = 50, min = 10, max = 200
                )
              ),
              shiny::actionButton("create_viz", "Create Visualization",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::selectInput(
                "img_format_results", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_viz", "Download Visualization",
                class = "btn-success"
              )
            ),
            shinydashboard::box(
              title = "Results Visualization", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::plotOutput("results_plot", height = "600px")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Enrichment Analysis
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "enrichment",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Enrichment Analysis Settings", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectInput("enrichment_contrast_select",
                "Select Contrast:",
                choices = NULL
              ),
              shiny::hr(),
              shiny::h4("Custom Enrichment Sets"),
              shiny::fileInput("custom_enrichment_file",
                "Upload Custom Sets CSV",
                accept = ".csv"
              ),
              shiny::helpText("CSV with columns: Lipid, Set_Name"),
              shiny::actionLink("info_custom_sets", " Help",
                icon = shiny::icon("info-circle")
              ),
              shiny::actionButton("load_custom_sets", "Load Custom Sets",
                class = "btn-primary"
              ),
              shiny::hr(),
              shiny::numericInput("min_set_size", "Min Set Size:", value = 5, min = 2),
              shiny::numericInput("max_set_size", "Max Set Size:", value = 500, min = 10),
              shiny::actionButton("run_enrichment", "Run Enrichment",
                class = "btn-primary"
              )
            ),
            shinydashboard::box(
              title = "Enrichment Results", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::selectInput("enrichment_contrast",
                "Display Contrast:",
                choices = NULL
              ),
              shiny::selectInput("enrichment_type", "Enrichment Type:",
                choices = c(
                  "Saturation" = "Saturation",
                  "Lipid Group" = "LipidGroup",
                  "Lipid Type" = "LipidType"
                )
              ),
              DT::dataTableOutput("enrichment_results_table"),
              shiny::downloadButton("download_enrichment",
                "Download Results",
                class = "btn-success"
              )
            )
          )
        ),

        # ------------------------------------------------------------------
        # Enrichment Visualization
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "enrichment_viz",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Enrichment Visualization Options", status = "primary",
              solidHeader = TRUE, width = 4,
              shiny::selectInput("enrichment_viz_contrast",
                "Select Contrast:",
                choices = NULL
              ),
              shiny::selectInput("enrichment_viz_type", "Enrichment Type:",
                choices = c(
                  "Saturation" = "Saturation",
                  "Lipid Group" = "LipidGroup",
                  "Lipid Type" = "LipidType"
                )
              ),
              shiny::selectInput("enrichment_plot_type", "Plot Type:",
                choices = c(
                  "Dot Plot" = "dotplot",
                  "Bar Plot" = "barplot"
                )
              ),
              shiny::numericInput("max_pathways", "Max Pathways to Show:",
                value = 15, min = 5, max = 30
              ),
              shiny::actionButton("create_enrichment_viz",
                "Create Visualization",
                class = "btn-primary"
              ),
              shiny::br(), shiny::br(),
              shiny::selectInput(
                "img_format_enrich", "Image Format:",
                c("PNG" = "png", "PDF" = "pdf")
              ),
              shiny::downloadButton("download_enrichment_viz",
                "Download Plot",
                class = "btn-success"
              )
            ),
            shinydashboard::box(
              title = "Enrichment Visualization", status = "info",
              solidHeader = TRUE, width = 8,
              shiny::plotOutput("enrichment_plot", height = "600px")
            )
          )
        ),

        # ------------------------------------------------------------------
        # Report Generation
        # ------------------------------------------------------------------
        shinydashboard::tabItem(
          tabName = "report",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Generate Analysis Report", status = "primary",
              solidHeader = TRUE, width = 6,
              shiny::textInput("report_title", "Report Title:",
                value = "Lipidomics Analysis Report"
              ),
              shiny::textInput("report_author", "Author:", value = ""),
              shiny::hr(),
              shiny::h4("Include Sections"),
              .checkbox_group_with_buttons(
                "report_sections", "Sections to include:",
                choices = c(
                  "Data Summary" = "data_summary",
                  "Normalization" = "normalization",
                  "Differential Analysis" = "diff_analysis",
                  "Enrichment Analysis" = "enrichment"
                ),
                selected = c(
                  "data_summary", "normalization",
                  "diff_analysis", "enrichment"
                )
              ),
              shiny::hr(),
              shiny::h4("Plot History"),
              shiny::verbatimTextOutput("plot_history_status"),
              shiny::hr(),
              shiny::h4("Output Format"),
              shiny::radioButtons("report_format", "Choose format:",
                choices  = c("HTML" = "html", "PDF" = "pdf"),
                selected = "html"
              ),
              shiny::helpText("Note: PDF requires LaTeX (tinytex). If PDF fails, use HTML."),
              shiny::br(),
              shiny::downloadButton("download_report",
                "Generate & Download Report",
                class = "btn-primary btn-lg"
              )
            ),
            shinydashboard::box(
              title = "Analysis Status", status = "info",
              solidHeader = TRUE, width = 6,
              shiny::verbatimTextOutput("report_status")
            )
          )
        )
      ) # end tabItems
    ) # end dashboardBody
  ) # end dashboardPage

  # ==========================================================================
  # Server
  # ==========================================================================
  server <- function(input, output, session) {
    options(shiny.sanitize.errors = TRUE)

    # ---- Reactive state ---------------------------------------------------
    values <- shiny::reactiveValues(
      raw_data = NULL,
      imputed_data = NULL,          # kept for legacy; use pre_impute_data
      pre_impute_data = NULL,       # normalized_data snapshot before imputation
      normalized_data = NULL,
      pre_batch_data = NULL,        # normalized_data snapshot before batch correction
      classification = NULL,
      custom_classification = NULL,
      use_custom_classification = FALSE,
      classification_data = NULL,
      diff_results = NULL,
      available_contrasts = NULL,
      custom_enrichment_sets = NULL,
      enrichment_results = NULL,

      # Single "current" plots (used for quick download)
      current_plot = NULL,
      current_raw_plot = NULL,
      current_norm_plot = NULL,
      current_results_plot = NULL,
      current_enrichment_plot = NULL,
      current_expression_plots = NULL,
      current_pipeline_plots = NULL,
      pipeline1_data = NULL,
      pipeline2_data = NULL,
      pipeline1_methods = NULL,
      pipeline2_methods = NULL,

      # Full plot history for the report
      plot_history = list()
    )

    # Convenience: show an error notification
    show_error <- function(e, context = "Operation") {
      shiny::showNotification(
        paste(context, "failed:", conditionMessage(e)),
        type = "error", duration = 10
      )
    }

    # Convenience: append a plot to the history list
    add_to_history <- function(section, label, plot_obj) {
      entry <- list(
        section = section, label = label, plot = plot_obj,
        time = Sys.time()
      )
      values$plot_history <- c(values$plot_history, list(entry))
    }

    # ---- Select All / Deselect All wiring --------------------------------
    .register_select_all_observers(
      session, input, "metadata_cols",
      function() c("Sample Name", "Sample Group", "Tumour ID", "Weight (mg)")
    )
    .register_select_all_observers(
      session, input, "norm_methods_1",
      function() get_normalization_methods()
    )
    .register_select_all_observers(
      session, input, "norm_methods_2",
      function() get_normalization_methods()
    )
    .register_select_all_observers(
      session, input, "groups_included",
      function() {
        shiny::req(values$raw_data, input$group_column)
        if (input$group_column %in% names(values$raw_data$metadata)) {
          unique(values$raw_data$metadata[[input$group_column]])
        } else {
          character(0)
        }
      }
    )
    .register_select_all_observers(
      session, input, "selected_samples",
      function() {
        shiny::req(values$raw_data)
        if ("Sample Name" %in% names(values$raw_data$metadata)) {
          values$raw_data$metadata$`Sample Name`
        } else {
          character(0)
        }
      }
    )
    .register_select_all_observers(
      session, input, "selected_expression_groups",
      function() {
        shiny::req(values$raw_data)
        gc <- if ("Sample Group" %in% names(values$raw_data$metadata)) {
          "Sample Group"
        } else {
          names(values$raw_data$metadata)[1]
        }
        if (!is.null(values$raw_data$metadata[[gc]])) {
          unique(values$raw_data$metadata[[gc]])
        } else {
          character(0)
        }
      }
    )
    .register_select_all_observers(
      session, input, "report_sections",
      function() c("data_summary", "normalization", "diff_analysis", "enrichment")
    )
    .register_select_all_observers(
      session, input, "raw_filter_groups",
      function() {
        shiny::req(values$raw_data)
        if ("Sample Group" %in% names(values$raw_data$metadata)) {
          unique(values$raw_data$metadata$`Sample Group`)
        } else {
          character(0)
        }
      }
    )
    .register_select_all_observers(
      session, input, "norm_filter_groups",
      function() {
        shiny::req(values$raw_data)
        if ("Sample Group" %in% names(values$raw_data$metadata)) {
          unique(values$raw_data$metadata$`Sample Group`)
        } else {
          character(0)
        }
      }
    )

    # ---- Load example data -----------------------------------------------
    shiny::observeEvent(input$load_example, {
      tryCatch(
        {
          example_data <- generate_example_data()
          values$raw_data <- load_lipidomics_data_from_df(
            example_data,
            input$metadata_cols
          )

          if (!is.null(values$raw_data$metadata$`Sample Name`)) {
            rownames(values$raw_data$numeric_data) <-
              values$raw_data$metadata$`Sample Name`
            rownames(values$raw_data$metadata) <-
              values$raw_data$metadata$`Sample Name`
          }

          lipid_names <- colnames(values$raw_data$numeric_data)
          values$classification <- classify_lipids(lipid_names)
          values$classification_data <- values$classification

          .update_ui_after_load(session, values)
          shiny::showNotification("Example data loaded successfully!")
        },
        error = function(e) show_error(e, "Loading example data")
      )
    })

    # ---- Load data from file ---------------------------------------------
    shiny::observeEvent(input$load_data, {
      shiny::req(input$file)
      tryCatch(
        {
          loaded <- load_lipidomics_data(
            input$file$datapath,
            input$metadata_cols
          )
          values$raw_data <- list(
            data         = loaded$data,
            metadata     = loaded$metadata,
            numeric_data = as.matrix(loaded$numeric_data)
          )

          if ("Sample Name" %in% names(values$raw_data$metadata)) {
            rownames(values$raw_data$numeric_data) <-
              values$raw_data$metadata$`Sample Name`
            rownames(values$raw_data$metadata) <-
              values$raw_data$metadata$`Sample Name`
          }

          lipid_names <- colnames(values$raw_data$numeric_data)
          values$classification <- classify_lipids(lipid_names)
          values$classification_data <- values$classification
          values$use_custom_classification <- FALSE

          .update_ui_after_load(session, values)
          shiny::showNotification("Data loaded successfully!", type = "message")
        },
        error = function(e) {
          shiny::showNotification(paste("Error loading data:", e$message),
            type = "error"
          )
        }
      )
    })

    # ---- Data summary / preview ------------------------------------------
    output$data_summary <- shiny::renderText({
      shiny::req(values$raw_data)
      paste0(
        "Samples: ", nrow(values$raw_data$data), "\n",
        "Features: ", ncol(values$raw_data$numeric_data), "\n",
        "Groups: ", paste(unique(values$raw_data$metadata$`Sample Group`),
          collapse = ", "
        )
      )
    })

    output$data_preview <- DT::renderDataTable({
      shiny::req(values$raw_data)
      DT::datatable(utils::head(values$raw_data$data),
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })

    # ---- Classification --------------------------------------------------
    output$classification_table <- DT::renderDataTable({
      cls <- if (values$use_custom_classification && !is.null(values$custom_classification)) {
        values$custom_classification
      } else {
        values$classification
      }
      if (!is.null(cls)) {
        DT::datatable(cls, options = list(scrollX = TRUE, pageLength = 10))
      } else {
        NULL
      }
    })

    output$download_classification <- shiny::downloadHandler(
      filename = function() .dl_name("lipid_classification", ext = "csv"),
      content = function(file) {
        cls <- if (values$use_custom_classification) {
          values$custom_classification
        } else {
          values$classification
        }
        utils::write.csv(cls, file, row.names = FALSE)
      }
    )

    shiny::observeEvent(input$load_custom_classification, {
      shiny::req(input$custom_classification_file)
      tryCatch(
        {
          values$custom_classification <-
            load_custom_classification(input$custom_classification_file$datapath)
          values$classification_data <- values$custom_classification
          values$use_custom_classification <- TRUE
          shiny::showNotification("Custom classification loaded!", type = "message")
        },
        error = function(e) show_error(e, "Loading custom classification")
      )
    })

    # Reset to automatic classification
    shiny::observeEvent(input$reset_classification, {
      shiny::req(values$classification)
      values$custom_classification <- NULL
      values$use_custom_classification <- FALSE
      values$classification_data <- values$classification
      # Restore default color_column choices
      shiny::updateSelectInput(session, "color_column",
        choices  = c("LipidGroup", "LipidType", "Saturation"),
        selected = "LipidGroup"
      )
      shiny::showNotification("Automatic classification restored.", type = "message")
    })

    # Dynamically update color_column choices when classification data changes
    shiny::observe({
      cls <- values$classification_data
      if (!is.null(cls)) {
        cls_cols <- setdiff(colnames(cls), "Lipid")
        if (length(cls_cols) > 0) {
          shiny::updateSelectInput(session, "color_column",
            choices  = cls_cols,
            selected = cls_cols[1]
          )
          # Also update enrichment type choices if custom columns present
          enrich_choices <- c()
          if ("Saturation" %in% cls_cols) enrich_choices <- c(enrich_choices, "Saturation" = "Saturation")
          if ("LipidGroup" %in% cls_cols) enrich_choices <- c(enrich_choices, "Lipid Group" = "LipidGroup")
          if ("LipidType" %in% cls_cols) enrich_choices <- c(enrich_choices, "Lipid Type" = "LipidType")
          # Add any custom columns not in the default set
          custom_cols <- setdiff(cls_cols, c("Saturation", "LipidGroup", "LipidType"))
          if (length(custom_cols) > 0) {
            extra <- custom_cols
            names(extra) <- custom_cols
            enrich_choices <- c(enrich_choices, extra)
          }
          if (length(enrich_choices) > 0) {
            for (id in c("enrichment_type", "enrichment_viz_type")) {
              shiny::updateSelectInput(session, id, choices = enrich_choices)
            }
          }
        }
      }
    })

    shiny::observeEvent(input$info_classification, {
      shiny::showModal(shiny::modalDialog(
        title = "Lipid Classification",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML("
          <p>Lipids are automatically classified based on their name prefixes:</p>
          <ul>
            <li><strong>LipidGroup</strong>: broad class (e.g., Glycerophospholipids, Sphingolipids)</li>
            <li><strong>LipidType</strong>: specific subclass (e.g., Phosphatidylcholine)</li>
            <li><strong>Saturation</strong>: SFA (0 double bonds), MUFA (1), PUFA (&gt;1)</li>
          </ul>
          <p>Upload a custom CSV (<em>Lipid</em> column + classification columns) to override this.</p>
        ")
      ))
    })

    # ---- Lipid / sample choice updates -----------------------------------
    shiny::observe({
      shiny::req(values$raw_data)
      if (!is.null(values$raw_data$numeric_data)) {
        shiny::updateSelectizeInput(session, "selected_lipids",
          choices = colnames(values$raw_data$numeric_data)
        )
      }
    })

    shiny::observe({
      shiny::req(values$raw_data)
      md <- values$raw_data$metadata
      if ("Sample Name" %in% names(md)) {
        shiny::updateCheckboxGroupInput(session, "selected_samples",
          choices  = md$`Sample Name`,
          selected = md$`Sample Name`
        )
      }
      gc <- if ("Sample Group" %in% names(md)) "Sample Group" else names(md)[1]
      if (!is.null(md[[gc]])) {
        grps <- unique(md[[gc]])
        shiny::updateCheckboxGroupInput(session, "selected_expression_groups",
          choices = grps, selected = grps
        )
      }
    })

    # ---- group_column -> groups_included update -------------------------
    shiny::observe({
      shiny::req(values$raw_data, input$group_column)
      md <- values$raw_data$metadata
      if (input$group_column %in% names(md)) {
        grps <- unique(md[[input$group_column]])
        grps <- grps[!is.na(grps)]
        shiny::updateCheckboxGroupInput(session, "groups_included",
          choices = grps, selected = grps
        )
      }
    })

    # ---- Raw data visualization -----------------------------------------
    shiny::observeEvent(input$create_raw_plot, {
      shiny::req(values$raw_data)
      tryCatch(
        {
          # Apply group filter
          raw_data_filtered <- values$raw_data
          if (!is.null(input$raw_filter_groups) && length(input$raw_filter_groups) > 0 &&
            "Sample Group" %in% names(values$raw_data$metadata)) {
            mask <- values$raw_data$metadata$`Sample Group` %in% input$raw_filter_groups
            raw_data_filtered <- list(
              numeric_data = values$raw_data$numeric_data[mask, , drop = FALSE],
              metadata     = values$raw_data$metadata[mask, , drop = FALSE]
            )
          }
          metadata <- if (isTRUE(input$raw_color_by_group)) raw_data_filtered$metadata else NULL

          if (input$plot_type == "heatmap") {
            hm <- create_heatmap_robust(
              data_matrix  = t(raw_data_filtered$numeric_data),
              metadata     = raw_data_filtered$metadata,
              group_column = "Sample Group",
              top_n        = input$raw_heatmap_top_n,
              title        = "Raw Data Heatmap - Top Variable Lipids"
            )
            values$current_raw_plot <- hm
            values$current_plot     <- hm
            add_to_history("Raw Data", "Raw heatmap", hm)
            output$raw_plot_ui <- shiny::renderUI(
              shiny::plotOutput("raw_heatmap_plot", height = "600px")
            )
            output$raw_heatmap_plot <- shiny::renderPlot({
              if (inherits(hm, "pheatmap")) grid::grid.draw(hm$gtable) else print(hm)
            })
          } else {
            top_n <- if (input$view_mode == "lipid") {
              input$top_n_lipids
            } else if (isTRUE(input$limit_top_samples)) {
              input$top_n_samples
            } else {
              NULL
            }
            plot  <- visualize_raw_data_improved(
              raw_data_filtered, input$plot_type, input$view_mode,
              top_n, metadata = metadata
            )
            values$current_raw_plot <- plot
            values$current_plot     <- plot
            add_to_history("Raw Data", paste("Raw", input$plot_type, input$view_mode), plot)
            output$raw_plot_ui <- shiny::renderUI(
              plotly::plotlyOutput("raw_plot", height = "500px")
            )
            output$raw_plot <- plotly::renderPlotly(plotly::ggplotly(plot))
          }
        },
        error = function(e) show_error(e, "Raw data visualization")
      )
    })

    # ---- Normalization pipeline comparison -------------------------------
    # Store pipeline data when button is clicked
    shiny::observeEvent(input$compare_pipelines, {
      shiny::req(values$raw_data)
      tryCatch(
        {
          values$pipeline1_data <- apply_normalizations(
            values$raw_data$numeric_data,
            input$norm_methods_1
          )
          values$pipeline2_data <- apply_normalizations(
            values$raw_data$numeric_data,
            input$norm_methods_2
          )
          values$pipeline1_methods <- input$norm_methods_1
          values$pipeline2_methods <- input$norm_methods_2
          shiny::showNotification("Pipeline comparison complete!", type = "message")
        },
        error = function(e) show_error(e, "Pipeline comparison")
      )
    })

    # Reactive renderPlot: re-builds automatically whenever data OR plot type changes
    output$pipeline_comparison <- shiny::renderPlot({
      shiny::req(values$pipeline1_data, values$pipeline2_data)
      md <- if (isTRUE(input$norm_color_by_group)) values$raw_data$metadata else NULL
      pt <- input$norm_compare_plot_type
      lbl1 <- if (!is.null(values$pipeline1_methods)) {
        paste(values$pipeline1_methods, collapse = " -> ")
      } else {
        "Pipeline 1"
      }
      lbl2 <- if (!is.null(values$pipeline2_methods)) {
        paste(values$pipeline2_methods, collapse = " -> ")
      } else {
        "Pipeline 2"
      }
      cmp_view <- if (is.null(input$norm_compare_view_mode)) {
        "sample"
      } else {
        input$norm_compare_view_mode
      }
      cmp_top_n <- if (isTRUE(input$norm_compare_limit_top)) {
        input$norm_compare_top_n
      } else {
        NULL
      }
      p1 <- create_pipeline_plot(values$pipeline1_data,
        title = paste("Pipeline 1:", lbl1),
        metadata = md, plot_type = pt,
        view_mode = cmp_view, top_n = cmp_top_n
      )
      p2 <- create_pipeline_plot(values$pipeline2_data,
        title = paste("Pipeline 2:", lbl2),
        metadata = md, plot_type = pt,
        view_mode = cmp_view, top_n = cmp_top_n
      )
      shiny::isolate(values$current_pipeline_plots <- list(p1 = p1, p2 = p2))
      gridExtra::grid.arrange(p1, p2, ncol = 1)
    })

    # ---- Apply normalization --------------------------------------------
    shiny::observeEvent(input$apply_normalization, {
      shiny::req(values$raw_data, input$chosen_pipeline)
      tryCatch(
        {
          methods <- if (input$chosen_pipeline == "1") {
            input$norm_methods_1
          } else {
            input$norm_methods_2
          }

          norm_mat <- apply_normalizations(values$raw_data$numeric_data, methods)

          if (!is.null(rownames(values$raw_data$numeric_data))) {
            rownames(norm_mat) <- rownames(values$raw_data$numeric_data)
          }

          values$normalized_data <- list(
            data         = cbind(values$raw_data$metadata, norm_mat),
            metadata     = values$raw_data$metadata,
            numeric_data = norm_mat
          )
          shiny::showNotification("Normalization applied successfully!")
        },
        error = function(e) show_error(e, "Applying normalization")
      )
    })

    # Normalization method descriptions modal
    shiny::observeEvent(input$info_norm_methods, {
      descs <- get_normalization_descriptions()
      html_rows <- paste(
        vapply(names(descs), function(m) {
          paste0("<dt><strong>", m, "</strong></dt><dd>", descs[[m]], "</dd>")
        }, character(1)),
        collapse = "\n"
      )
      shiny::showModal(shiny::modalDialog(
        title = "Normalization Method Descriptions",
        size = "l",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML(paste0(
          "<dl class='dl-horizontal'>", html_rows, "</dl>",
          "<hr><p><strong>Note on Quantile Normalization:</strong> ",
          "After quantile normalization, per-sample boxplots will look nearly ",
          "identical - this is the <em>expected, correct</em> behaviour, not a bug. ",
          "The method forces every sample to share the same intensity distribution.</p>"
        ))
      ))
    })

    # Ellipse help modal
    shiny::observeEvent(input$info_ellipses, {
      shiny::showModal(shiny::modalDialog(
        title = "Ellipse Types",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML("
          <dl>
            <dt><strong>None</strong></dt>
            <dd>No ellipse drawn. Only the scatter points are shown.</dd>
            <dt><strong>Confidence Ellipse</strong></dt>
            <dd>A normal-distribution confidence region at the chosen level (default 95%).
                Useful to assess whether groups are statistically separable.</dd>
            <dt><strong>Visual Circle</strong></dt>
            <dd>A t-distribution ellipse at 80% level. Gives a rough visual grouping
                without strong statistical assumptions.</dd>
          </dl>
          <p><em>Tip:</em> Sample labels and group ellipses are independent.
          Enable 'Show Sample Labels' to see individual sample names in addition to
          the group ellipses.</p>
        ")
      ))
    })

    # ---- Normalized data visualization ----------------------------------
    shiny::observeEvent(input$create_norm_plot, {
      shiny::req(values$normalized_data)
      tryCatch(
        {
          group_col <- input$group_column
          if (is.null(rownames(values$normalized_data$metadata))) {
            rownames(values$normalized_data$metadata) <-
              rownames(values$normalized_data$numeric_data)
          }

          # Filter by selected groups
          # Filter by norm_filter_groups (independent of PCA groups_included)
          filter_grps <- if (length(input$norm_filter_groups) > 0) {
            input$norm_filter_groups
          } else {
            NULL
          }
          if (!is.null(filter_grps) && group_col %in% colnames(values$normalized_data$metadata)) {
            mask <- values$normalized_data$metadata[[group_col]] %in% filter_grps
            filt_dat <- values$normalized_data$numeric_data[mask, , drop = FALSE]
            filt_md <- values$normalized_data$metadata[mask, , drop = FALSE]
          } else if (!is.null(input$groups_included) && length(input$groups_included) > 0 &&
            group_col %in% colnames(values$normalized_data$metadata)) {
            mask <- values$normalized_data$metadata[[group_col]] %in% input$groups_included
            filt_dat <- values$normalized_data$numeric_data[mask, , drop = FALSE]
            filt_md <- values$normalized_data$metadata[mask, , drop = FALSE]
          } else {
            filt_dat <- values$normalized_data$numeric_data
            filt_md <- values$normalized_data$metadata
          }

          if (input$norm_plot_type %in% c("boxplot", "violin", "density")) {
            dl <- list(numeric_data = filt_dat, metadata = filt_md)
            md_col <- if (isTRUE(input$norm_color_by_group_viz)) filt_md else NULL
            norm_top_n <- if (isTRUE(input$norm_limit_top_samples)) {
              input$norm_top_n_samples
            } else {
              NULL
            }
            plot <- visualize_raw_data_improved(
              dl, input$norm_plot_type, "sample",
              top_n = norm_top_n,
              metadata = md_col, group_column = group_col
            )
            values$current_norm_plot <- plot
            add_to_history("Normalized Data", paste("Norm", input$norm_plot_type), plot)
            output$norm_plot_ui <- shiny::renderUI(
              plotly::plotlyOutput("norm_plot", height = "500px")
            )
            output$norm_plot <- plotly::renderPlotly(plotly::ggplotly(plot))

          } else if (input$norm_plot_type == "heatmap") {
            hm <- create_heatmap_robust(
              data_matrix  = t(filt_dat),
              metadata     = filt_md,
              group_column = group_col,
              top_n        = input$norm_heatmap_top_n,
              title        = "Normalised Data Heatmap - Top Variable Lipids"
            )
            values$current_norm_plot <- hm
            add_to_history("Normalized Data", "Norm heatmap", hm)
            output$norm_plot_ui <- shiny::renderUI(
              shiny::plotOutput("norm_heatmap_plot", height = "600px")
            )
            output$norm_heatmap_plot <- shiny::renderPlot({
              if (inherits(hm, "pheatmap")) grid::grid.draw(hm$gtable) else print(hm)
            })

          } else if (input$norm_plot_type == "pca") {
            shiny::req(group_col)
            pca_res <- perform_pca(filt_dat, filt_md, group_col)
            plot <- create_pca_plot_with_ellipses(
              pca_data = pca_res$pca_data,
              variance_explained = pca_res$variance_explained,
              ellipse_type = input$ellipse_type,
              show_sample_labels = isTRUE(input$show_sample_labels)
            )
            values$current_norm_plot <- plot
            add_to_history("Normalized Data", "Norm PCA", plot)
            output$norm_plot_ui <- shiny::renderUI(
              plotly::plotlyOutput("norm_plot", height = "500px")
            )
            output$norm_plot <- plotly::renderPlotly(plotly::ggplotly(plot))

          } else if (input$norm_plot_type == "plsda") {
            shiny::req(group_col)
            plsda_res <- perform_plsda(as.matrix(filt_dat), filt_md, group_col)
            plot <- create_plsda_plot_with_ellipses(
              plsda_data = plsda_res$scores_data,
              ellipse_type = input$ellipse_type,
              show_sample_labels = isTRUE(input$show_sample_labels)
            )
            values$current_norm_plot <- plot
            add_to_history("Normalized Data", "Norm PLS-DA", plot)
            output$norm_plot_ui <- shiny::renderUI(
              plotly::plotlyOutput("norm_plot", height = "500px")
            )
            output$norm_plot <- plotly::renderPlotly(plotly::ggplotly(plot))
          }
        },
        error = function(e) show_error(e, "Creating normalized data plot")
      )
    })

    # ---- Lipid expression -----------------------------------------------
    shiny::observeEvent(input$create_expression_plot, {
      shiny::req(values$raw_data, input$selected_lipids)
      tryCatch(
        {
          if (input$expression_data_type == "normalized" &&
            !is.null(values$normalized_data)) {
            data_use <- values$normalized_data$numeric_data
            md_use <- values$normalized_data$metadata
          } else {
            data_use <- values$raw_data$numeric_data
            md_use <- values$raw_data$metadata
          }

          sel_samples <- if (input$expression_selection_mode == "samples") {
            input$selected_samples
          } else {
            NULL
          }
          sel_groups <- if (input$expression_selection_mode == "groups") {
            input$selected_expression_groups
          } else {
            NULL
          }

          gc <- if ("Sample Group" %in% names(md_use)) {
            "Sample Group"
          } else {
            names(md_use)[1]
          }

          plots <- create_lipid_expression_barplot(
            data_matrix      = data_use,
            metadata         = md_use,
            selected_lipids  = input$selected_lipids,
            selected_samples = sel_samples,
            selected_groups  = sel_groups,
            group_column     = gc,
            data_type        = input$expression_data_type
          )

          values$current_expression_plots <- plots

          # Add each lipid plot to history
          if (is.list(plots) && !inherits(plots, "ggplot")) {
            for (nm in names(plots)) {
              add_to_history("Lipid Expression", paste("Expression", nm), plots[[nm]])
            }
          } else {
            add_to_history("Lipid Expression", "Expression Plot", plots)
          }

          shiny::showNotification("Plots created!", type = "message")
        },
        error = function(e) show_error(e, "Creating expression plot")
      )
    })

    output$expression_plots_ui <- shiny::renderUI({
      shiny::req(values$current_expression_plots)
      if (is.list(values$current_expression_plots) &&
        !inherits(values$current_expression_plots, "ggplot")) {
        do.call(shiny::tagList, lapply(
          names(values$current_expression_plots),
          function(lipid) {
            shinydashboard::box(
              title = lipid, width = 12,
              shiny::plotOutput(
                paste0("expr_plot_", make.names(lipid)),
                height = "400px"
              )
            )
          }
        ))
      } else {
        shiny::plotOutput("expr_plot_single", height = "400px")
      }
    })

    shiny::observe({
      plots <- values$current_expression_plots
      if (is.list(plots) && !inherits(plots, "ggplot")) {
        lapply(names(plots), function(lipid) {
          local({
            my_lipid <- lipid
            output_name <- paste0("expr_plot_", make.names(my_lipid))
            output[[output_name]] <- shiny::renderPlot(plots[[my_lipid]])
          })
        })
      } else if (!is.null(plots)) {
        output$expr_plot_single <- shiny::renderPlot(plots)
      }
    })

    # ---- Differential analysis contrasts --------------------------------
    shiny::observe({
      md <- if (!is.null(values$normalized_data)) {
        values$normalized_data$metadata
      } else if (!is.null(values$raw_data)) {
        values$raw_data$metadata
      } else {
        NULL
      }
      shiny::req(md, input$group_col_diff)
      if (input$group_col_diff %in% names(md)) {
        grps <- levels(factor(md[[input$group_col_diff]]))
        values$available_contrasts <- create_default_contrasts(grps)
      }
    })

    output$contrast_selection_ui <- shiny::renderUI({
      shiny::req(values$available_contrasts)
      .checkbox_group_with_buttons(
        "selected_contrasts", "Choose Contrasts:",
        choices = values$available_contrasts,
        selected = values$available_contrasts
      )
    })

    # Wire Select/Deselect for dynamically rendered contrast checkboxes
    shiny::observeEvent(input$selected_contrasts_select_all, {
      shiny::updateCheckboxGroupInput(session, "selected_contrasts",
        selected = values$available_contrasts
      )
    })
    shiny::observeEvent(input$selected_contrasts_deselect_all, {
      shiny::updateCheckboxGroupInput(session, "selected_contrasts",
        selected = character(0)
      )
    })

    # ---- Info modals for diff analysis -----------------------------------
    shiny::observeEvent(input$info_methods, {
      shiny::showModal(shiny::modalDialog(
        title = "Differential Analysis Methods", easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML("
          <h4>limma</h4>
          <p>Linear Models for Microarray and Omics Data. Fast, robust for
             continuous normalized data. Recommended for lipidomics.</p>
          <h4>EdgeR</h4>
          <p>Originally for RNA-seq count data; uses negative-binomial models.
             Data are rounded to integer counts before fitting.</p>
        ")
      ))
    })

    shiny::observeEvent(input$info_contrasts, {
      shiny::showModal(shiny::modalDialog(
        title = "Contrast Selection", easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML("
          <h4>Automatic pairwise contrasts</h4>
          <p>Each contrast compares two groups: <code>GroupA - GroupB</code>
          tests whether GroupA differs from GroupB. Positive logFC means
          higher in GroupA. All pairwise contrasts are generated automatically.</p>
          <h4>Custom contrasts</h4>
          <p>Use the text box to define complex comparisons in standard
          limma notation. Group names must exactly match those in the
          selected group column (use <code>make.names()</code> rules:
          spaces become dots). One contrast per line.</p>
          <p><strong>Examples for a 2x2 design (Diet x Treatment):</strong></p>
          <pre style='font-size:12px'>
# Diet effect (collapsed across treatment):
(ND_Vehicle + ND_DabTram)/2 - (HFD_Vehicle + HFD_DabTram)/2

# Treatment effect (collapsed across diet):
(ND_DabTram + HFD_DabTram)/2 - (ND_Vehicle + HFD_Vehicle)/2

# Interaction (does treatment effect differ by diet?):
(HFD_DabTram - HFD_Vehicle) - (ND_DabTram - ND_Vehicle)</pre>
          <p><em>Tip:</em> if group names contain spaces or special
          characters, they are internally sanitised with
          <code>make.names()</code> - use the sanitised form in custom
          contrasts (e.g. <code>A375.ND.DabTram</code> for
          <code>A375 ND DabTram</code>).</p>
        ")
      ))
    })

    # ---- Run differential analysis ---------------------------------------
    shiny::observeEvent(input$run_diff_analysis, {
      # At least one of: checkbox contrasts or custom text contrasts must be present
      checkbox_contrasts <- input$selected_contrasts
      custom_raw <- trimws(input$custom_contrasts_text)
      custom_contrasts <- if (nzchar(custom_raw)) {
        lines <- trimws(strsplit(custom_raw, "\n")[[1]])
        lines[nzchar(lines) & !startsWith(lines, "#")]
      } else {
        character(0)
      }
      all_contrasts <- unique(c(checkbox_contrasts, custom_contrasts))

      if (length(all_contrasts) == 0) {
        shiny::showNotification(
          "Select at least one contrast or enter a custom contrast.",
          type = "warning", duration = 6
        )
        return()
      }

      shiny::req(values$normalized_data)
      tryCatch(
        {
          md <- values$normalized_data$metadata

          results <- perform_differential_analysis(
            data_matrix    = values$normalized_data$numeric_data,
            metadata       = md,
            group_column   = input$group_col_diff,
            contrasts_list = all_contrasts,
            method         = input$diff_method
          )
          values$diff_results <- results

          contrast_choices <- names(results$results)
          for (id in c(
            "contrast_display", "contrast_select",
            "enrichment_contrast_select",
            "enrichment_contrast", "enrichment_viz_contrast"
          )) {
            shiny::updateSelectInput(session, id,
              choices  = contrast_choices,
              selected = contrast_choices[1]
            )
          }
          shiny::showNotification(
            paste("Differential analysis completed using", input$diff_method),
            type = "message"
          )
        },
        error = function(e) show_error(e, "Differential analysis")
      )
    })

    output$diff_summary <- shiny::renderText({
      if (is.null(values$diff_results)) {
        return("No results yet. Run the analysis first.")
      }

      out <- paste0(
        "Method: ", values$diff_results$method, "\n",
        "Contrasts: ", length(values$diff_results$results), "\n\n"
      )
      for (nm in names(values$diff_results$results)) {
        res <- values$diff_results$results[[nm]]
        n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
        n_up <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
        n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
        out <- paste0(
          out,
          "=== ", nm, " ===\n",
          "  Total: ", nrow(res), "\n",
          "  Significant (adj.P < 0.05): ", n_sig, "\n",
          "    Increased abundance: ", n_up, "\n",
          "    Decreased abundance: ", n_down, "\n\n"
        )
      }
      out
    })

    output$diff_results_table <- DT::renderDataTable(
      {
        shiny::req(values$diff_results, input$contrast_display)
        tryCatch(
          {
            res <- values$diff_results$results[[input$contrast_display]]
            if (is.null(res)) {
              return(DT::datatable(data.frame(Message = "No results for this contrast")))
            }
            if (!is.data.frame(res)) res <- as.data.frame(res)
            if (!"Lipid" %in% colnames(res)) {
              res <- data.frame(Lipid = rownames(res), res, stringsAsFactors = FALSE)
            }

            num_cols <- vapply(res, is.numeric, logical(1))
            res[num_cols] <- lapply(res[num_cols], round, 4)

            DT::datatable(res,
              options = list(
                pageLength = 15, scrollX = TRUE,
                scrollY = "400px", searching = TRUE
              ),
              rownames = FALSE,
              filter = "top",
              selection = "multiple"
            )
          },
          error = function(e) {
            DT::datatable(data.frame(Error = paste("Display error:", e$message)))
          }
        )
      },
      server = TRUE
    )

    # ---- Results visualization ------------------------------------------
    shiny::observeEvent(input$create_viz, {
      shiny::req(values$diff_results, input$contrast_select)
      tryCatch(
        {
          sel_res <- values$diff_results$results[[input$contrast_select]]

          if (input$viz_type == "volcano") {
            cls <- if (isTRUE(input$color_by_class)) values$classification_data else NULL
            clr_by <- if (isTRUE(input$color_by_class)) input$color_column else NULL

            plot <- create_volcano_plot_labeled(
              sel_res,
              title = paste("Volcano Plot:", input$contrast_select),
              logfc_threshold = input$logfc_threshold,
              pval_threshold = input$pval_threshold,
              top_labels = input$top_labels,
              classification_data = cls,
              color_by = clr_by
            )
            values$current_results_plot <- plot
            add_to_history(
              "Results",
              paste("Volcano", input$contrast_select), plot
            )
            output$results_plot <- shiny::renderPlot(plot)
          } else if (input$viz_type == "heatmap") {
            sig_idx <- which(sel_res$adj.P.Val < 0.05 & abs(sel_res$logFC) > 1)

            if (length(sig_idx) == 0) {
              output$results_plot <- shiny::renderPlot(
                ggplot2::ggplot() +
                  ggplot2::annotate("text",
                    x = 0.5, y = 0.5,
                    label = "No significant features found",
                    size = 6
                  ) +
                  ggplot2::theme_void()
              )
            } else {
              n_feat <- min(length(sig_idx), input$heatmap_top_n)
              feat_names <- rownames(sel_res)[sig_idx[seq_len(n_feat)]]
              avail <- intersect(
                feat_names,
                colnames(values$normalized_data$numeric_data)
              )

              if (length(avail) == 0) {
                output$results_plot <- shiny::renderPlot(
                  ggplot2::ggplot() +
                    ggplot2::annotate("text",
                      x = 0.5, y = 0.5,
                      label = paste0(
                        "No significant lipid names matched the data matrix.\n",
                        "Check that differential analysis was run on the\n",
                        "same normalised data currently loaded."
                      ),
                      size = 5
                    ) +
                    ggplot2::theme_void()
                )
              } else {
                # Safely filter to the contrast groups.
                # Only attempt simple "A - B" contrasts; for complex expressions
                # (custom contrasts with parentheses/arithmetic) show all samples.
                group_vals <- as.character(
                  values$normalized_data$metadata[[input$group_col_diff]]
                )
                is_simple <- !grepl("[()*/+]", input$contrast_select)
                mask <- if (is_simple) {
                  parts <- trimws(strsplit(input$contrast_select, " - ")[[1]])
                  group_vals %in% parts
                } else {
                  rep(TRUE, length(group_vals))
                }
                # Require at least 2 samples after filtering; otherwise show all
                if (sum(mask) < 2L) mask <- rep(TRUE, length(group_vals))

                filt_dat <- values$normalized_data$numeric_data[mask, avail, drop = FALSE]
                filt_md  <- values$normalized_data$metadata[mask, , drop = FALSE]

                hm <- create_heatmap_robust(
                  t(filt_dat), filt_md, input$group_col_diff,
                  top_n = length(avail),
                  title = paste("Heatmap:", input$contrast_select)
                )
                values$current_results_plot <- hm
                output$results_plot <- shiny::renderPlot({
                  if (inherits(hm, "pheatmap")) {
                    grid::grid.draw(hm$gtable)
                  } else {
                    hm
                  }
                })
              }
            }
          }
        },
        error = function(e) show_error(e, "Creating visualization")
      )
    })

    # ---- Enrichment ------------------------------------------------------
    shiny::observeEvent(input$load_custom_sets, {
      shiny::req(input$custom_enrichment_file)
      tryCatch(
        {
          values$custom_enrichment_sets <-
            load_custom_enrichment_sets(input$custom_enrichment_file$datapath)
          shiny::showNotification("Custom enrichment sets loaded!", type = "message")
        },
        error = function(e) show_error(e, "Loading custom enrichment sets")
      )
    })

    shiny::observeEvent(input$info_custom_sets, {
      shiny::showModal(shiny::modalDialog(
        title = "Custom Enrichment Sets", easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::HTML("
          <p>Upload a CSV with two columns:</p>
          <ul>
            <li><code>Lipid</code>: must match the lipid names in your data exactly</li>
            <li><code>Set_Name</code>: name of the functional set/pathway</li>
          </ul>
          <p>A lipid can appear in multiple rows to belong to several sets.</p>
        ")
      ))
    })

    shiny::observeEvent(input$run_enrichment, {
      shiny::req(values$diff_results, values$classification)
      tryCatch(
        {
          cls <- if (values$use_custom_classification) {
            values$custom_classification
          } else {
            values$classification
          }

          values$enrichment_results <- perform_enrichment_analysis(
            results_list = values$diff_results$results,
            classification_data = cls,
            min_set_size = input$min_set_size,
            max_set_size = input$max_set_size,
            custom_sets = values$custom_enrichment_sets
          )
          shiny::showNotification("Enrichment analysis completed!", type = "message")
        },
        error = function(e) {
          if (grepl("corrupt|fgsea", e$message, ignore.case = TRUE)) {
            shiny::showNotification(
              paste(
                "fgsea error - try reinstalling: install.packages('fgsea')\n",
                e$message
              ),
              type = "error", duration = 15
            )
          } else {
            show_error(e, "Enrichment analysis")
          }
        }
      )
    })

    output$enrichment_results_table <- DT::renderDataTable({
      shiny::req(
        values$enrichment_results,
        input$enrichment_contrast, input$enrichment_type
      )
      tryCatch(
        {
          cr <- values$enrichment_results[[input$enrichment_contrast]]
          if (!is.null(cr) && input$enrichment_type %in% names(cr)) {
            ed <- cr[[input$enrichment_type]]
            if (nrow(ed) > 0) {
              DT::datatable(ed, options = list(scrollX = TRUE))
            } else {
              DT::datatable(data.frame(Message = "No results for this type"))
            }
          } else {
            DT::datatable(data.frame(Message = "Enrichment type not available"))
          }
        },
        error = function(e) {
          DT::datatable(data.frame(Error = paste("Display error:", e$message)))
        }
      )
    })

    shiny::observeEvent(input$create_enrichment_viz, {
      shiny::req(
        values$enrichment_results,
        input$enrichment_viz_contrast, input$enrichment_viz_type
      )
      tryCatch(
        {
          cr <- values$enrichment_results[[input$enrichment_viz_contrast]]
          if (!is.null(cr) && input$enrichment_viz_type %in% names(cr)) {
            ed <- cr[[input$enrichment_viz_type]]

            plot <- if (nrow(ed) == 0) {
              ggplot2::ggplot() +
                ggplot2::annotate("text",
                  x = 0.5, y = 0.5,
                  label = "No enrichment results"
                ) +
                ggplot2::theme_void()
            } else if (input$enrichment_plot_type == "dotplot") {
              create_enrichment_dotplot(
                ed,
                title = paste(
                  "Enrichment:", input$enrichment_viz_contrast,
                  "-", input$enrichment_viz_type
                ),
                max_pathways = input$max_pathways
              )
            } else {
              create_enrichment_barplot(
                ed,
                title = paste(
                  "Enrichment:", input$enrichment_viz_contrast,
                  "-", input$enrichment_viz_type
                ),
                max_pathways = input$max_pathways
              )
            }

            values$current_enrichment_plot <- plot
            add_to_history(
              "Enrichment",
              paste(
                "Enrichment", input$enrichment_viz_contrast,
                input$enrichment_viz_type
              ), plot
            )
            output$enrichment_plot <- shiny::renderPlot(plot)
          }
        },
        error = function(e) show_error(e, "Enrichment visualization")
      )
    })

    # ---- Preprocessing: missing value summary ----------------------------
    output$missing_value_summary <- shiny::renderText({
      if (is.null(values$normalized_data)) {
        return("Apply normalisation first, then return here to impute.")
      }
      m       <- values$normalized_data$numeric_data
      n_miss  <- sum(is.na(m))
      pct     <- round(100 * n_miss / length(m), 2)
      imputed <- !is.null(values$pre_impute_data)
      paste0(
        "Total values  : ", length(m), "\n",
        "Missing (NA)  : ", n_miss, " (", pct, "%)\n",
        "Imputation    : ", if (imputed) "Applied" else "Not yet applied"
      )
    })

    # ---- Preprocessing: imputation ---------------------------------------
    shiny::observeEvent(input$run_imputation, {
      if (is.null(values$normalized_data)) {
        shiny::showNotification(
          "Please apply normalisation first, then return here to impute.",
          type = "warning", duration = 6
        )
        return()
      }
      tryCatch({
        k_val <- if (input$imputation_method == "knn") input$imputation_k else 5L

        # Snapshot for reset
        values$pre_impute_data <- values$normalized_data

        imp <- impute_missing_values(
          values$normalized_data$numeric_data,
          method = input$imputation_method,
          k      = k_val
        )

        md <- values$normalized_data$metadata
        values$normalized_data <- list(
          data         = cbind(md, imp),
          metadata     = md,
          numeric_data = imp
        )
        shiny::showNotification(
          paste("Imputation applied to normalised data:", input$imputation_method),
          type = "message"
        )
      }, error = function(e) show_error(e, "Imputation"))
    })

    shiny::observeEvent(input$reset_imputation, {
      if (is.null(values$pre_impute_data)) {
        shiny::showNotification("No imputation to reset.", type = "warning")
        return()
      }
      values$normalized_data  <- values$pre_impute_data
      values$pre_impute_data  <- NULL
      shiny::showNotification(
        "Imputation reset. Pre-imputation normalised data restored.",
        type = "message"
      )
    })

    shiny::observeEvent(input$info_imputation, {
      descs    <- get_imputation_descriptions()
      html_rows <- paste(
        vapply(names(descs), function(m) {
          paste0("<dt><strong>", m, "</strong></dt><dd>", descs[[m]], "</dd>")
        }, character(1)),
        collapse = "\n"
      )
      shiny::showModal(shiny::modalDialog(
        title     = "Imputation Method Descriptions",
        size      = "l",
        easyClose = TRUE,
        footer    = shiny::modalButton("Close"),
        shiny::HTML(paste0("<dl class='dl-horizontal'>", html_rows, "</dl>"))
      ))
    })

    # ---- Preprocessing: batch correction --------------------------------
    shiny::observeEvent(input$run_batch_correction, {
      # -- Guard 1: normalisation must be applied first ----------------------
      if (is.null(values$normalized_data)) {
        shiny::showNotification(
          "Please apply normalisation first (Normalisation tab) before batch correction.",
          type = "warning", duration = 6
        )
        return()
      }

      # -- Guard 2: batch column must exist in metadata ----------------------
      md <- values$normalized_data$metadata
      if (is.null(input$batch_column) || !nzchar(input$batch_column) ||
          !input$batch_column %in% colnames(md)) {
        shiny::showNotification(
          paste0("No valid batch column selected. Add a batch column to your ",
                 "CSV (e.g. 'Batch', 'Run', 'Plate') and reload the data."),
          type = "warning", duration = 8
        )
        return()
      }

      # -- Guard 3: batch column must have 2...(n-1) unique values -----------
      batch_vals <- md[[input$batch_column]]
      n_batches  <- length(unique(batch_vals[!is.na(batch_vals)]))
      n_samples  <- nrow(md)

      if (n_batches < 2L) {
        shiny::showNotification(
          paste0("The selected batch column ('", input$batch_column,
                 "') has only one unique value. ",
                 "Choose a column where samples are labelled by their batch ",
                 "(e.g. run date, plate ID, instrument run)."),
          type = "error", duration = 10
        )
        return()
      }
      if (n_batches == n_samples) {
        shiny::showNotification(
          paste0("Every sample has a unique value in '", input$batch_column,
                 "'. Batch correction requires multiple samples per batch. ",
                 "Select a column that groups samples into 2 or more shared batches."),
          type = "error", duration = 10
        )
        return()
      }

      # -- Guard 4: warn if batch column is same as group column -------------
      # This is statistically questionable (removes biology alongside batch),
      # but we allow it with a warning. The collinearity check inside
      # correct_batch_effects() will handle complete confounding gracefully.
      if (!is.null(input$batch_group_column) &&
          input$batch_column == input$batch_group_column) {
        shiny::showNotification(
          paste0("Note: batch column and group column are the same ('",
                 input$batch_column, "'). Group protection will be disabled ",
                 "to avoid collinearity. Consider adding a dedicated batch ",
                 "column (e.g. run date or plate ID) to your data file."),
          type = "warning", duration = 10
        )
      }

      tryCatch({
        # Snapshot metadata before any assignment (avoids reactive mid-read)
        snap_md  <- values$normalized_data$metadata
        snap_mat <- values$normalized_data$numeric_data

        # Store full snapshot for reset
        values$pre_batch_data <- values$normalized_data

        corrected_mat <- correct_batch_effects(
          data_matrix  = snap_mat,
          metadata     = snap_md,
          batch_column = input$batch_column,
          group_column = input$batch_group_column,
          method       = input$batch_method
        )

        values$normalized_data <- list(
          data         = cbind(snap_md, corrected_mat),
          metadata     = snap_md,
          numeric_data = corrected_mat
        )

        # correct_batch_effects() silently falls back from "combat" to "limma"
        # when the sva package isn't installed (with only a console warning());
        # reflect the method actually used, not the one requested.
        method_used <- attr(corrected_mat, "method_used")
        if (!is.null(method_used) && method_used != input$batch_method) {
          shiny::showNotification(
            paste0(
              "ComBat requires the 'sva' package, which is not installed. ",
              "Batch correction was applied using limma::removeBatchEffect instead. ",
              "Install 'sva' (BiocManager::install(\"sva\")) to use ComBat."
            ),
            type = "warning", duration = 15
          )
        } else {
          shiny::showNotification(
            paste0("Batch correction applied using ", input$batch_method,
                   ". Verify the result with a PCA plot."),
            type = "message", duration = 6
          )
        }
      }, error = function(e) {
        # Roll back snapshot so normalized_data stays valid
        values$normalized_data <- values$pre_batch_data
        values$pre_batch_data  <- NULL
        shiny::showNotification(
          paste0("Batch correction failed: ", conditionMessage(e),
                 ". Check that batch and group columns are not perfectly confounded."),
          type = "error", duration = 12
        )
      })
    })

    shiny::observeEvent(input$reset_batch_correction, {
      shiny::req(values$pre_batch_data)
      values$normalized_data <- values$pre_batch_data
      values$pre_batch_data  <- NULL
      shiny::showNotification(
        "Batch correction reset. Pre-batch normalised data restored.",
        type = "message"
      )
    })

    shiny::observeEvent(input$info_batch, {
      shiny::showModal(shiny::modalDialog(
        title     = "Batch Effect Correction",
        easyClose = TRUE,
        footer    = shiny::modalButton("Close"),
        shiny::HTML("
          <h4>limma - removeBatchEffect</h4>
          <p>Fits a linear model that includes batch as a covariate and
          returns residuals without the batch component. Fast and reliable
          for most designs. Requires the batch structure to be fully known
          and balanced.</p>
          <h4>ComBat (sva package)</h4>
          <p>Parametric empirical Bayes adjustment. Generally more robust
          when batch effects are large or when batch sizes are unequal.
          Requires the <code>sva</code> Bioconductor package.</p>
          <hr/>
          <p><strong>Important:</strong> always apply batch correction
          <em>after</em> normalisation and verify the result with a PCA
          plot (Normalised Data Visualization tab). Samples should cluster by
          biological group, not by batch, after correction.</p>
        ")
      ))
    })

    # ---- Preprocessing status -------------------------------------------
    output$preprocessing_status <- shiny::renderText({
      paste0(
        "Data loaded      : ", !is.null(values$raw_data), "\n",
        "Normalisation    : ",
        if (!is.null(values$normalized_data)) "Applied" else "Not applied", "\n",
        "Imputation       : ",
        if (!is.null(values$pre_impute_data)) "Applied" else "Not applied", "\n",
        "Batch correction : ",
        if (!is.null(values$pre_batch_data)) "Applied" else "Not applied"
      )
    })

    # ---- Report status ---------------------------------------------------
    output$report_status <- shiny::renderText({
      paste0(
        "Data loaded:             ", !is.null(values$raw_data), "\n",
        "Normalization applied:   ", !is.null(values$normalized_data), "\n",
        "Diff. analysis done:     ", !is.null(values$diff_results), "\n",
        "Enrichment done:         ", !is.null(values$enrichment_results), "\n",
        "Plots in history:        ", length(values$plot_history)
      )
    })

    # Plot history summary shown on report tab
    output$plot_history_status <- shiny::renderText({
      n <- length(values$plot_history)
      if (n == 0) {
        return("No plots generated yet.")
      }
      entries <- vapply(values$plot_history, function(e) {
        paste0("[", e$section, "] ", e$label)
      }, character(1))
      paste(c(paste(n, "plot(s) will be included:"), entries), collapse = "\n")
    })

    # ==========================================================================
    # Download handlers
    # ==========================================================================

    # Descriptive filenames for all downloads

    output$download_raw_plot <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_raw, "pdf")) "pdf" else "png"
        .dl_name("raw", input$view_mode, input$plot_type, ext = ext)
      },
      content = function(file) {
        shiny::req(values$current_raw_plot)
        .save_plot(values$current_raw_plot, file, input$img_format_raw)
      }
    )

    output$download_norm_plot <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_norm, "pdf")) "pdf" else "png"
        .dl_name("normalized", input$norm_plot_type,
          if (!is.null(input$ellipse_type)) input$ellipse_type else "",
          ext = ext
        )
      },
      content = function(file) {
        shiny::req(values$current_norm_plot)
        .save_plot(values$current_norm_plot, file, input$img_format_norm)
      }
    )

    output$download_pipeline_comparison <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_pipeline, "pdf")) "pdf" else "png"
        .dl_name("pipeline_comparison",
          paste(input$norm_methods_1, collapse = "-"), "vs",
          paste(input$norm_methods_2, collapse = "-"),
          ext = ext
        )
      },
      content = function(file) {
        shiny::req(values$raw_data)
        tryCatch(
          {
            dl_view <- if (is.null(input$norm_compare_view_mode)) {
              "sample"
            } else {
              input$norm_compare_view_mode
            }
            dl_top_n <- if (isTRUE(input$norm_compare_limit_top)) {
              input$norm_compare_top_n
            } else {
              NULL
            }
            p1 <- create_pipeline_plot(
              apply_normalizations(values$raw_data$numeric_data, input$norm_methods_1),
              paste("Pipeline 1:", paste(input$norm_methods_1, collapse = "->")),
              metadata = if (input$norm_color_by_group) values$raw_data$metadata else NULL,
              plot_type = input$norm_compare_plot_type,
              view_mode = dl_view, top_n = dl_top_n
            )
            p2 <- create_pipeline_plot(
              apply_normalizations(values$raw_data$numeric_data, input$norm_methods_2),
              paste("Pipeline 2:", paste(input$norm_methods_2, collapse = "->")),
              metadata = if (input$norm_color_by_group) values$raw_data$metadata else NULL,
              plot_type = input$norm_compare_plot_type,
              view_mode = dl_view, top_n = dl_top_n
            )
            combined <- gridExtra::arrangeGrob(p1, p2, ncol = 1)
            ggplot2::ggsave(file, combined,
              width = 12, height = 10, units = "in",
              device = if (identical(input$img_format_pipeline, "pdf")) {
                "pdf"
              } else {
                "png"
              },
              dpi = 300
            )
          },
          error = function(e) show_error(e, "Pipeline comparison download")
        )
      }
    )

    output$download_viz <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_results, "pdf")) "pdf" else "png"
        .dl_name("results", input$viz_type,
          gsub("[^A-Za-z0-9]", "_", input$contrast_select),
          ext = ext
        )
      },
      content = function(file) {
        shiny::req(values$current_results_plot)
        tryCatch(
          {
            is_pdf <- identical(input$img_format_results, "pdf")
            if (inherits(values$current_results_plot, "pheatmap")) {
              if (is_pdf) {
                grDevices::pdf(file,
                  width = 12, height = 10,
                  useDingbats = FALSE
                )
              } else {
                grDevices::png(file,
                  width = 12, height = 10,
                  units = "in", res = 300
                )
              }
              grid::grid.draw(values$current_results_plot$gtable)
              grDevices::dev.off()
            } else {
              ggplot2::ggsave(file, values$current_results_plot,
                width = 12, height = 10, units = "in",
                device = if (is_pdf) "pdf" else "png", dpi = 300
              )
            }
          },
          error = function(e) show_error(e, "Results download")
        )
      }
    )

    output$download_current_table <- shiny::downloadHandler(
      filename = function() {
        .dl_name("diff", gsub("[^A-Za-z0-9]", "_", input$contrast_display),
          ext = "csv"
        )
      },
      content = function(file) {
        shiny::req(values$diff_results, input$contrast_display)
        tryCatch(
          {
            res <- values$diff_results$results[[input$contrast_display]]
            if (!"Lipid" %in% colnames(res)) {
              res <- data.frame(Lipid = rownames(res), res, stringsAsFactors = FALSE)
            }
            utils::write.csv(res, file, row.names = FALSE)
          },
          error = function(e) show_error(e, "Table download")
        )
      }
    )

    # Excel export: truncate sheet names to Excel's 31-char limit
    output$download_results <- shiny::downloadHandler(
      filename = function() .dl_name("differential_all", ext = "xlsx"),
      content = function(file) {
        shiny::req(values$diff_results)
        tryCatch(
          {
            wb <- openxlsx::createWorkbook()
            for (nm in names(values$diff_results$results)) {
              sheet_nm <- .truncate_sheet_name(nm) # Excel sheet names are capped at 31 chars
              openxlsx::addWorksheet(wb, sheet_nm)
              res <- values$diff_results$results[[nm]]
              if (!"Lipid" %in% colnames(res)) {
                res <- data.frame(Lipid = rownames(res), res, stringsAsFactors = FALSE)
              }
              openxlsx::writeData(wb, sheet_nm, res, rowNames = FALSE)
            }
            openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
          },
          error = function(e) show_error(e, "Excel download")
        )
      }
    )

    output$download_enrichment <- shiny::downloadHandler(
      filename = function() {
        .dl_name("enrichment",
          gsub("[^A-Za-z0-9]", "_", input$enrichment_contrast),
          input$enrichment_type,
          ext = "csv"
        )
      },
      content = function(file) {
        shiny::req(
          values$enrichment_results,
          input$enrichment_contrast, input$enrichment_type
        )
        tryCatch(
          {
            ed <- values$enrichment_results[[input$enrichment_contrast]][[
              input$enrichment_type
            ]]
            utils::write.csv(ed, file, row.names = FALSE)
          },
          error = function(e) show_error(e, "Enrichment download")
        )
      }
    )

    output$download_enrichment_viz <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_enrich, "pdf")) "pdf" else "png"
        .dl_name("enrichment_plot",
          input$enrichment_plot_type,
          gsub("[^A-Za-z0-9]", "_", input$enrichment_viz_contrast),
          input$enrichment_viz_type,
          ext = ext
        )
      },
      content = function(file) {
        shiny::req(values$current_enrichment_plot)
        .save_plot(values$current_enrichment_plot, file, input$img_format_enrich)
      }
    )

    output$download_expression_plots <- shiny::downloadHandler(
      filename = function() {
        ext <- if (identical(input$img_format_expression, "pdf")) "pdf" else "png"
        lipid_tag <- if (!is.null(input$selected_lipids) &&
          length(input$selected_lipids) > 0) {
          paste(make.names(input$selected_lipids[seq_len(min(
            3,
            length(input$selected_lipids)
          ))]), collapse = "-")
        } else {
          "no_lipid"
        }
        .dl_name("expression", lipid_tag, ext = ext)
      },
      content = function(file) {
        shiny::req(values$current_expression_plots)
        tryCatch(
          {
            is_pdf <- identical(input$img_format_expression, "pdf")
            if (is.list(values$current_expression_plots) &&
              !inherits(values$current_expression_plots, "ggplot")) {
              combined <- gridExtra::arrangeGrob(
                grobs = values$current_expression_plots, ncol = 1
              )
              ggplot2::ggsave(file, combined,
                width = 12,
                height = 6 * length(values$current_expression_plots),
                units = "in",
                device = if (is_pdf) "pdf" else "png", dpi = 300
              )
            } else {
              .save_plot(
                values$current_expression_plots, file,
                input$img_format_expression
              )
            }
          },
          error = function(e) show_error(e, "Expression plot download")
        )
      }
    )

    # ---- Report generation -----------------------------------------------
    output$download_report <- shiny::downloadHandler(
      filename = function() {
        ext <- if (input$report_format == "pdf") "pdf" else "html"
        paste0("lipidomics_report_", Sys.Date(), ".", ext)
      },
      content = function(file) {
        prog <- shiny::showNotification(
          "Generating report... please wait.",
          duration = NULL, type = "message"
        )
        on.exit(try(shiny::removeNotification(prog), silent = TRUE))

        tryCatch(
          {
            report_dir <- file.path(
              tempdir(),
              paste0("rpt_", format(Sys.time(), "%Y%m%d%H%M%S"))
            )
            dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
            temp_rmd <- file.path(report_dir, "report.Rmd")

            # --- Save current plots as named files (primary route) ---
            # This works even if the user never clicked the explicit buttons;
            # any plot that was ever created and stored is included.
            .save_rpt_plot <- function(plot_obj, filename) {
              fp <- file.path(report_dir, filename)
              tryCatch(
                {
                  if (inherits(plot_obj, "pheatmap")) {
                    grDevices::png(fp, width = 10, height = 8, units = "in", res = 150)
                    grid::grid.draw(plot_obj$gtable)
                    grDevices::dev.off()
                  } else if (inherits(plot_obj, "ggplot")) {
                    ggplot2::ggsave(fp, plot_obj,
                      width = 10, height = 6,
                      dpi = 150, units = "in"
                    )
                  } else {
                    return(NULL)
                  }
                  filename # return filename on success
                },
                error = function(e) {
                  if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
                  NULL
                }
              )
            }

            # Build simple named list of plot file paths
            pf <- list()
            if (!is.null(values$current_raw_plot)) {
              pf$raw_plot <- .save_rpt_plot(values$current_raw_plot, "raw_data.png")
            }
            if (!is.null(values$current_pipeline_plots)) {
              combined <- gridExtra::arrangeGrob(
                values$current_pipeline_plots$p1,
                values$current_pipeline_plots$p2,
                ncol = 1
              )
              fp <- file.path(report_dir, "pipeline_comparison.png")
              tryCatch(
                ggplot2::ggsave(fp, combined,
                  width = 10, height = 10,
                  dpi = 150, units = "in"
                ),
                error = function(e) NULL
              )
              pf$pipeline_plot <- "pipeline_comparison.png"
            }
            if (!is.null(values$current_norm_plot)) {
              pf$norm_plot <- .save_rpt_plot(values$current_norm_plot, "normalized.png")
            }
            if (!is.null(values$current_results_plot)) {
              pf$results_plot <- .save_rpt_plot(values$current_results_plot, "results.png")
            }
            if (!is.null(values$current_enrichment_plot)) {
              pf$enrichment_plot <- .save_rpt_plot(values$current_enrichment_plot, "enrichment.png")
            }

            # Extra plots from session history (appended as appendix)
            extra_plots <- list()
            for (entry in values$plot_history) {
              fn <- paste0(make.names(paste0(entry$section, "_", entry$label)), ".png")
              saved <- .save_rpt_plot(entry$plot, fn)
              if (!is.null(saved)) {
                extra_plots[[length(extra_plots) + 1]] <- list(
                  file = fn, label = entry$label, section = entry$section
                )
              }
            }

            rmd <- .build_report_rmd_with_plots(
              title              = input$report_title,
              author             = input$report_author,
              sections           = input$report_sections,
              raw_data           = values$raw_data,
              normalized_data    = values$normalized_data,
              diff_results       = values$diff_results,
              enrichment_results = values$enrichment_results,
              output_format      = input$report_format,
              plot_files         = pf,
              extra_plots        = extra_plots
            )

            # Open, write, and explicitly close the connection: writeLines()
            # does not close connections passed in already-open, and without
            # an explicit close() the write buffer may not be flushed to disk
            # before rmarkdown::render() reads the file below, truncating it
            # mid-content.
            rmd_con <- file(temp_rmd, "w", encoding = "UTF-8")
            writeLines(rmd, con = rmd_con)
            close(rmd_con)

            fmt <- if (input$report_format == "pdf") {
              rmarkdown::pdf_document(toc = TRUE)
            } else {
              rmarkdown::html_document(toc = TRUE, toc_float = TRUE, theme = "flatly")
            }

            out_file <- tryCatch(
              {
                rmarkdown::render(temp_rmd,
                  output_format = fmt,
                  output_dir = report_dir, quiet = TRUE
                )
              },
              error = function(pdf_err) {
                if (input$report_format == "pdf") {
                  shiny::showNotification(
                    "PDF failed - generating HTML instead. Install tinytex for PDF support.",
                    type = "warning", duration = 10
                  )
                  rmarkdown::render(temp_rmd,
                    output_format = rmarkdown::html_document(
                      toc = TRUE, toc_float = TRUE, theme = "flatly"
                    ),
                    output_dir = report_dir, quiet = TRUE
                  )
                } else {
                  stop(pdf_err)
                }
              }
            )

            file.copy(out_file, file, overwrite = TRUE)
            shiny::showNotification("Report generated successfully!", type = "message")
          },
          error = function(e) {
            shiny::showNotification(paste("Report failed:", e$message),
              type = "error", duration = 15
            )
            writeLines(paste("# Report Error\n\n", e$message), file)
          }
        )
      }
    )
  } # end server

  shiny::shinyApp(
    ui = ui, server = server,
    options = list(port = port)
  )
}


# ===========================================================================
# Internal helpers
# ===========================================================================

# Save a ggplot or pheatmap to a file (PNG or PDF)
.save_plot <- function(plot_obj, file, fmt) {
  is_pdf <- identical(fmt, "pdf")
  if (inherits(plot_obj, "pheatmap")) {
    if (is_pdf) {
      grDevices::pdf(file, width = 12, height = 10, useDingbats = FALSE)
    } else {
      grDevices::png(file, width = 12, height = 10, units = "in", res = 300)
    }
    on.exit(if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE),
      add = TRUE
    )
    grid::grid.draw(plot_obj$gtable)
  } else if (inherits(plot_obj, "ggplot")) {
    ggplot2::ggsave(file, plot_obj,
      width = 12, height = 8, units = "in",
      device = if (is_pdf) "pdf" else "png", dpi = 300
    )
  } else if (inherits(plot_obj, c("grob", "gtable"))) {
    if (is_pdf) {
      grDevices::pdf(file, width = 12, height = 10, useDingbats = FALSE)
    } else {
      grDevices::png(file, width = 12, height = 10, units = "in", res = 300)
    }
    on.exit(if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE),
      add = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(plot_obj)
  } else {
    stop("Cannot save object of class: ", paste(class(plot_obj), collapse = ", "))
  }
}

# Update UI elements after data load (shared between file upload and example)
.update_ui_after_load <- function(session, values) {
  md <- values$raw_data$metadata
  cols <- names(md)

  shiny::updateSelectInput(session, "group_column",
    choices = cols,
    selected = if ("Sample Group" %in% cols) {
      "Sample Group"
    } else {
      cols[1]
    }
  )
  shiny::updateSelectInput(session, "group_col_diff",
    choices = cols,
    selected = if ("Sample Group" %in% cols) {
      "Sample Group"
    } else {
      cols[1]
    }
  )
  # Batch correction selectors
  # "Sample Name" is never a valid batch variable (every sample is unique)
  # so exclude it from choices. Also exclude numeric-ish columns like Weight.
  batch_candidates <- cols[!cols %in% c("Sample Name")]
  batch_default <- if (length(batch_candidates) > 0) {
    # Prefer a column that is neither Sample Name nor Sample Group
    non_bio <- batch_candidates[!batch_candidates %in% c("Sample Group")]
    if (length(non_bio) > 0) non_bio[1] else batch_candidates[1]
  } else {
    NULL
  }
  shiny::updateSelectInput(session, "batch_column",
    choices  = batch_candidates,
    selected = batch_default
  )
  shiny::updateSelectInput(session, "batch_group_column",
    choices  = cols,
    selected = if ("Sample Group" %in% cols) "Sample Group" else cols[1]
  )

  if ("Sample Name" %in% cols) {
    sn <- md$`Sample Name`
    shiny::updateCheckboxGroupInput(session, "selected_samples",
      choices = sn, selected = sn
    )
  }
  if ("Sample Group" %in% cols) {
    grps <- unique(md$`Sample Group`)
    shiny::updateCheckboxGroupInput(session, "selected_expression_groups",
      choices = grps, selected = grps
    )
    # Update raw and norm viz group filters
    shiny::updateCheckboxGroupInput(session, "raw_filter_groups",
      choices = grps, selected = grps
    )
    shiny::updateCheckboxGroupInput(session, "norm_filter_groups",
      choices = grps, selected = grps
    )
    shiny::updateCheckboxGroupInput(session, "groups_included",
      choices = grps, selected = grps
    )
  }
}


# ===========================================================================
# Report builder
# ===========================================================================

#' Build the Rmd Content for the Analysis Report
#'
#' Generates a complete R Markdown document as a character string.
#' \code{plot_files} is a simple named list of PNG file paths
#' (raw_plot, norm_plot, pipeline_plot, results_plot, enrichment_plot).
#' \code{extra_plots} is an optional list of additional plot entries
#' from the session history.
#'
#' @param title          Report title.
#' @param author         Author name.
#' @param sections       Character vector of section keys to include.
#' @param raw_data       Raw data list.
#' @param normalized_data Normalized data list.
#' @param diff_results   Differential analysis results list.
#' @param enrichment_results Enrichment analysis results list.
#' @param output_format  Either \code{"html"} or \code{"pdf"}.
#' @param plot_files     Named list of plot PNG paths
#'   (raw_plot, norm_plot, pipeline_plot, results_plot, enrichment_plot).
#' @param extra_plots    Optional list of additional plot entries from
#'   session history; each entry has \code{$file}, \code{$section},
#'   \code{$label}.
#' @return A single character string containing the complete Rmd document.
#' @noRd
.build_report_rmd_with_plots <- function(title, author, sections,
                                        raw_data, normalized_data,
                                        diff_results, enrichment_results,
                                        output_format = "html",
                                        plot_files = list(),
                                        extra_plots = list()) {
  # ---- YAML header -------------------------------------------------------
  yaml <- if (output_format == "pdf") {
    paste0(
      '---\ntitle: "', title, '"\nauthor: "', author, '"\n',
      'date: "', format(Sys.Date(), "%B %d, %Y"), '"\n',
      "output:\n  pdf_document:\n    toc: true\n    toc_depth: 2\n---\n\n"
    )
  } else {
    paste0(
      '---\ntitle: "', title, '"\nauthor: "', author, '"\n',
      'date: "', format(Sys.Date(), "%B %d, %Y"), '"\n',
      "output:\n  html_document:\n",
      "    toc: true\n    toc_float: true\n",
      "    theme: flatly\n    highlight: tango\n---\n\n"
    )
  }

  setup <- paste0(
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(\n",
    "  echo    = FALSE,\n",
    "  warning = FALSE,\n",
    "  message = FALSE,\n",
    "  fig.width  = 10,\n",
    "  fig.height = 6,\n",
    '  fig.align  = "center",\n',
    '  out.width  = "100%"\n',
    ")\n```\n\n"
  )

  # Helper: embed a plot image in the Rmd.
  # path is just the filename (PNG is in the same dir as the Rmd).
  # file.exists is NOT checked here because knitr resolves relative to Rmd dir.
  embed_plot <- function(path, caption = "", width = "100%") {
    if (!is.null(path) && nchar(path) > 0) {
      paste0(
        '\n```{r echo=FALSE, out.width="', width,
        '", fig.cap="', gsub('"', "'", caption), '"}\n',
        'knitr::include_graphics("', path, '")\n',
        "```\n\n"
      )
    } else {
      ""
    }
  }

  content <- ""

  # ---- 1. Data Summary ---------------------------------------------------
  if ("data_summary" %in% sections) {
    content <- paste0(content, "# Data Summary\n\n")

    if (!is.null(raw_data)) {
      n_s <- nrow(raw_data$data)
      n_f <- ncol(raw_data$numeric_data)
      content <- paste0(
        content,
        "## Dataset Overview\n\n",
        "| Metric | Value |\n|--------|-------|\n",
        "| Samples | ", n_s, " |\n",
        "| Lipid features | ", n_f, " |\n\n"
      )

      if ("Sample Group" %in% names(raw_data$metadata)) {
        gt <- table(raw_data$metadata$`Sample Group`)
        content <- paste0(
          content,
          "## Sample Distribution\n\n",
          "| Group | Count |\n|-------|-------|\n",
          paste(
            vapply(names(gt), function(g) {
              paste0("| ", g, " | ", gt[[g]], " |")
            }, character(1)),
            collapse = "\n"
          ),
          "\n\n"
        )
      }

      content <- paste0(
        content,
        embed_plot(plot_files$raw_plot, "Raw Data Distribution")
      )
    } else {
      content <- paste0(content, "_No data loaded._\n\n")
    }
  }

  # ---- 2. Normalisation --------------------------------------------------
  if ("normalization" %in% sections) {
    content <- paste0(content, "\n---\n\n# Normalisation\n\n")

    if (!is.null(normalized_data)) {
      content <- paste0(
        content,
        "Normalisation has been applied to the dataset.\n\n",
        embed_plot(plot_files$pipeline_plot, "Pipeline Comparison"),
        embed_plot(plot_files$norm_plot, "Normalised Data Distribution")
      )
    } else {
      content <- paste0(content, "_Normalisation has not been applied._\n\n")
    }
  }

  # ---- 3. Differential Analysis ------------------------------------------
  if ("diff_analysis" %in% sections) {
    content <- paste0(content, "\n---\n\n# Differential Analysis\n\n")

    if (!is.null(diff_results)) {
      content <- paste0(
        content,
        "- **Method:** ", diff_results$method, "\n",
        "- **Contrasts tested:** ", length(diff_results$results), "\n\n",
        "## Results Summary\n\n"
      )

      for (nm in names(diff_results$results)) {
        res <- diff_results$results[[nm]]
        n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
        n_up <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
        n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
        content <- paste0(
          content,
          "### ", nm, "\n\n",
          "| Metric | Count |\n|--------|-------|\n",
          "| Total features | ", nrow(res), " |\n",
          "| Significant (adj.P < 0.05) | ", n_sig, " |\n",
          "| Increased abundance | ", n_up, " |\n",
          "| Decreased abundance | ", n_down, " |\n\n"
        )
      }

      content <- paste0(
        content,
        embed_plot(plot_files$results_plot, "Differential Analysis Results")
      )

      # Top 10 features from first contrast
      if (length(diff_results$results) > 0) {
        first_nm <- names(diff_results$results)[1]
        first_res <- diff_results$results[[first_nm]]
        sig <- first_res[!is.na(first_res$adj.P.Val) &
          first_res$adj.P.Val < 0.05, ]
        if (nrow(sig) > 0) {
          sig <- sig[order(sig$adj.P.Val), ]
          top10 <- utils::head(sig, 10)
          content <- paste0(
            content,
            "## Top Features -- ", first_nm, "\n\n",
            "| Lipid | LogFC | adj.P.Val | Direction |\n",
            "|-------|-------|-----------|----------|\n",
            paste(vapply(seq_len(nrow(top10)), function(i) {
              dir <- if (top10$logFC[i] > 0) "Increased" else "Decreased"
              paste0(
                "| ", rownames(top10)[i],
                " | ", round(top10$logFC[i], 3),
                " | ", format(top10$adj.P.Val[i], scientific = TRUE, digits = 3),
                " | ", dir, " |"
              )
            }, character(1)), collapse = "\n"),
            "\n\n"
          )
        }
      }
    } else {
      content <- paste0(content, "_Differential analysis has not been run._\n\n")
    }
  }

  # ---- 4. Enrichment -----------------------------------------------------
  if ("enrichment" %in% sections) {
    content <- paste0(content, "\n---\n\n# Enrichment Analysis\n\n")

    if (!is.null(enrichment_results)) {
      content <- paste0(
        content,
        "**Contrasts analysed:** ", length(enrichment_results), "\n\n",
        embed_plot(plot_files$enrichment_plot, "Enrichment Analysis")
      )

      for (nm in names(enrichment_results)) {
        cr <- enrichment_results[[nm]]
        content <- paste0(content, "### ", nm, "\n\n")
        for (et in names(cr)) {
          ed <- cr[[et]]
          sig <- if (!is.null(ed) && nrow(ed) > 0) {
            ed[!is.na(ed$padj) & ed$padj < 0.25, ]
          } else {
            data.frame()
          }
          if (nrow(sig) > 0) {
            top <- utils::head(sig[order(sig$pval), ], 5)
            content <- paste0(
              content,
              "**", et, "** - significant pathways (padj < 0.25):\n\n",
              "| Pathway | NES | p-value | padj |\n",
              "|---------|-----|---------|------|\n",
              paste(vapply(
                seq_len(nrow(top)), function(i) {
                  paste0(
                    "| ", top$pathway[i],
                    " | ", round(top$NES[i], 3),
                    " | ", format(top$pval[i], scientific = TRUE, digits = 2),
                    " | ", format(top$padj[i], scientific = TRUE, digits = 2), " |"
                  )
                },
                character(1)
              ), collapse = "\n"),
              "\n\n"
            )
          }
        }
      }
    } else {
      content <- paste0(content, "_Enrichment analysis has not been run._\n\n")
    }
  }

  # ---- 5. Extra plots from session history --------------------------------
  if (length(extra_plots) > 0) {
    content <- paste0(content, "\n---\n\n# Additional Visualisations\n\n")
    for (e in extra_plots) {
      content <- paste0(
        content,
        "## ", e$label, "\n\n",
        embed_plot(e$file, e$label)
      )
    }
  }

  # ---- 6. Session info ---------------------------------------------------
  content <- paste0(
    content,
    "\n---\n\n# Session Information\n\n",
    "- **Report generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "- **R version:** ", R.version.string, "\n\n"
  )

  paste0(yaml, setup, content)
}
