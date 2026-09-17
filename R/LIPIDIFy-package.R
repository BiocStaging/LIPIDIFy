#' LIPIDIFy: Comprehensive Lipidomics Data Analysis with Interactive Visualization
#'
#' LIPIDIFy provides a comprehensive toolkit for end-to-end lipidomics data
#' analysis, including missing value imputation, batch effect correction,
#' normalization, differential abundance analysis using limma and edgeR,
#' gene set enrichment analysis, and extensive visualization capabilities.
#' Lipid names are automatically classified by class, subclass, and
#' fatty-acid saturation. The package offers both an interactive Shiny
#' interface for bench biologists (see \code{\link{launch_lipidomics_app}})
#' and fully scriptable R functions for bioinformaticians. It supports
#' flexible custom lipid classification schemes and user-defined enrichment
#' sets.
#'
#' @section Getting started:
#' See \code{vignette("LIPIDIFy")} for a full walkthrough of the
#' programmatic workflow.
#'
#' @keywords internal
"_PACKAGE"

## Every package listed under Imports in DESCRIPTION is also declared here
## with @importFrom, so that NAMESPACE records exactly which functions LIPIDIFy
## relies on. Code in R/ still calls them as pkg::fun() for clarity.
## Keep this list in sync when adding or removing a pkg::fun() call.
#' @importFrom dplyr everything
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom FactoMineR PCA
#' @importFrom fgsea fgseaMultilevel
#' @importFrom ggplot2 aes annotate element_text geom_bar geom_boxplot
#' @importFrom ggplot2 geom_density geom_histogram geom_hline geom_point
#' @importFrom ggplot2 geom_text geom_violin geom_vline ggplot ggplotGrob ggsave
#' @importFrom ggplot2 guide_legend guides labs scale_color_gradient
#' @importFrom ggplot2 scale_color_manual scale_fill_brewer scale_fill_manual
#' @importFrom ggplot2 scale_size_continuous stat_ellipse theme theme_minimal
#' @importFrom ggplot2 theme_void
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices dev.cur dev.off pdf png
#' @importFrom grid grid.draw grid.newpage
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom limma contrasts.fit eBayes lmFit makeContrasts removeBatchEffect
#' @importFrom limma topTable
#' @importFrom openxlsx addWorksheet createWorkbook saveWorkbook writeData
#' @importFrom pheatmap pheatmap
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom pls plsr
#' @importFrom rmarkdown html_document pdf_document render
#' @importFrom scales hue_pal
#' @importFrom shiny actionButton actionLink br checkboxGroupInput checkboxInput
#' @importFrom shiny column conditionalPanel downloadButton downloadHandler
#' @importFrom shiny fileInput fluidRow h4 h5 helpText hr HTML icon isolate
#' @importFrom shiny modalButton modalDialog numericInput observe observeEvent p
#' @importFrom shiny plotOutput radioButtons reactiveValues removeNotification
#' @importFrom shiny renderPlot renderText renderUI req selectInput
#' @importFrom shiny selectizeInput shinyApp showModal showNotification strong
#' @importFrom shiny tagList tags textAreaInput textInput uiOutput
#' @importFrom shiny updateCheckboxGroupInput updateSelectInput
#' @importFrom shiny updateSelectizeInput verbatimTextOutput
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardPage
#' @importFrom shinydashboard dashboardSidebar menuItem sidebarMenu tabItem
#' @importFrom shinydashboard tabItems
#' @importFrom stats median model.matrix reorder rlnorm rnorm setNames var
#' @importFrom stringr regex str_detect str_extract str_extract_all str_match
#' @importFrom stringr str_remove_all str_trim
#' @importFrom tidyr all_of pivot_longer
#' @importFrom utils capture.output globalVariables head read.csv write.csv
#' @importFrom withr with_preserve_seed with_seed
NULL
