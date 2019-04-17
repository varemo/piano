#' Explore GSA results
#'
#' Explore GSA results interactively in a web browser using \code{shiny}.
#'
#' Additional gene-level information, e.g. alternative names or description,
#' can be supplied via the \code{geneAnnot} argument. This information will
#' show up in the gene table and the gene summary tabs.
#'
#' @param gsaRes an object of class \code{GSAres}, as returned from
#' \code{runGSA()} or an object returned from \code{runGSAhyper()}.
#' @param browser a logical, whether or not to open the Shiny app in a browser
#' window. Set to \code{FALSE} to open an interactive window directly in
#' RStudio.
#' @param geneAnnot a \code{data.frame}, containing gene annotation. The first
#' column should be gene IDs matching those in \code{gsares}.
#' @param genesets a character vector or list (named or un-named) of character vectors
#' containing subsets of gene-set names that can be selected and displayed in the
#' network plot.
#' @return Does not return any object.
#' @author Leif Varemo \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{GSAheatmap}}, \code{\link{networkPlot2}}
#' @examples
#'
#'    # Load example input data to GSA:
#'    data("gsa_input")
#'
#'    # Load gene set collection:
#'    gsc <- loadGSC(gsa_input$gsc)
#'
#'    # Run gene set analysis:
#'    gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions,
#'                     gsc=gsc, nPerm=500)
#'
#'    # Explore results:
#'    \dontrun{exploreGSAres(gsares)}
#'
exploreGSAres <- function(gsaRes, browser=TRUE, geneAnnot=NULL, genesets) {

  gsares <- gsaRes

  # Argument checking:
  if(!is(gsares, "GSAres")) stop("argument gsares is not of class GSAres")
  if(!is.logical(browser)) stop("argument browser is not a logical")
  if(!is.null(geneAnnot)) {
    if(sum(duplicated(geneAnnot[,1]))>0) stop("geneAnnot may not contain duplicated gene IDs")
    if(sum(geneAnnot[,1]%in%rownames(gsares$geneLevelStats)) == 0) stop("no overlap between genes in geneAnnot and gsares")
  }
  if(missing(genesets)) genesets <- list("No available gene-sets" = "No available gene-sets")
  if(is(genesets, "list")) {
    if(unique(unlist(lapply(genesets, class)))[1] != "character" | length(unique(unlist(lapply(genesets, class))))>1) stop("genesets has to be a character vector or a list of character vectors")
  } else {
    if(!is(genesets, "character")) stop("genesets has to be a character vector or a list of character vectors")
    genesets <- list(genesets)
  }
  if(is.null(names(genesets))) names(genesets) <- paste("Gene-set list", seq(from=1, to=length(genesets)))
  if(any(names(genesets) == "")) names(genesets)[names(genesets)==""] <- paste("Gene-set list", seq(from=1, to=sum(names(genesets)=="")))

  # Check packages now, otherwise delayed error during app browsing, if missing:
  if(!"shiny"%in%installed.packages()[,"Package"]) stop("missing package shiny")
  if(!"shinyjs"%in%installed.packages()[,"Package"]) stop("missing package shinyjs")
  if(!"DT"%in%installed.packages()[,"Package"]) stop("missing package DT")

  # Set needed variables for network plot:
  selectable_classes_network <- list()
  if(!all(is.na(gsares$pDistinctDirUp)) | !all(is.na(gsares$pDistinctDirDn))) {
     selectable_classes_network$"Distinct directional" <- "distinct_both"
  }
  if(!all(is.na(gsares$pNonDirectional))) {
    selectable_classes_network$"Non directional" <- "non_NULL"
  }
  if(!all(is.na(gsares$pMixedDirUp)) | !all(is.na(gsares$pMixedDirDn))) {
    selectable_classes_network$"Mixed directional up" <- "mixed_up"
    selectable_classes_network$"Mixed directional down" <- "mixed_down"
  }

  selected_class_network <- selectable_classes_network[[1]]

  if(selected_class_network == "distinct_both") {
    cutoff_network <- signif(sort(c(gsares$pAdjDistinctDirUp,gsares$pAdjDistinctDirDn), decreasing=FALSE)[20],2)
  } else if(selected_class_network == "non_NULL") {
    cutoff_network <- signif(sort(gsares$pAdjNonDirectional, decreasing=FALSE)[20],2)
  } else if(selected_class_network == "mixed_up") {
    cutoff_network <- signif(sort(c(gsares$pAdjMixedDirUp), decreasing=FALSE)[20],2)
  }
  if(cutoff_network==0) cutoff_network <- 1e-6


  # App object with ui and server
  app <- shinyApp(

# ===================================================================================
    # ui

#     ui <- fluidPage(
ui <- dashboardPage(
      dashboardHeader(disable=TRUE),
      dashboardSidebar(disable=TRUE),
      dashboardBody(navbarPage("Explore your GSA results",id="navbarpage",selected="tab_gsainfo", collapsible=TRUE,

      # fluidRow(
      #   column(8,
      #          h2(HTML("Explore gene-set analysis results"))
      #          )#,
      #   #column(4,
      #   #       div(imageOutput(system.file("shiny", "loggo.png", package="piano"), height="70px"), style="text-align: right;")
      #   #)
      # ),
      #HTML("<br>"),
      #fluidRow(
        #column(12,
          #tabsetPanel(id="tabset1",
            #type = "tabs",
            tabPanel(title="Run info", value="tab_gsainfo",
                     useShinyjs(),
                     #includeCSS("css/cosmo.css"),
                     #includeCSS("css/toggle_switch.css"),
                     #includeCSS("css/style.css"),
                     includeCSS(system.file("shiny", "css", "cosmo.css", package="piano")),
                     includeCSS(system.file("shiny", "css", "toggle_switch.css", package="piano")),
                     includeCSS(system.file("shiny", "css", "style.css", package="piano")),
                     tags$head(
                       tags$style(
                         HTML(".shiny-notification {
                    position:fixed;
                    top: calc(25%);
                    left: calc(25%);
                    width: 300px;
                    }"
                         )
                       )
                     ),
                     div(id="loading_tab_gsainfo",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                         ),
                     hidden(div(id="tab_gsainfo",
                     HTML("<br>"),
                     column(5,
                            box(title="Gene/gene-set info", status="primary", width="100%",
                                      htmlOutput("text_runinfo1")
                            ),
                            box(title="Analysis details", status="primary", width="100%",
                                      htmlOutput("text_runinfo2")
                            ),
                            wellPanel(style="background-color: #E2F3FF;",
                                      htmlOutput("text_runinfo3")
                            )
                     ),
                     column(7,
                            tabBox(width="100%",
                                   tabPanel("Gene-set size", value="gene_set_size_hist",
                                            plotOutput("gsc_boxplot1", height="340px", width="100%"),
                                            HTML("<div style='max-width: 600px'><b>Number of <font color='#009900'>genes</font> per <font color='#9966ff'>gene-set</font>.</b>
                          Gives you an idea about sizes of the gene-sets used in the analysis.
                          Are you looking at small specific gene-sets with few genes, or large
                          general gene-sets? Note that these stats are <em>after</em> filtering the
                          input GSC according to available gene-level statistics and the gsSizeLim parameter.
                          <br><br></div>")
                                   ),
                                   tabPanel("Gene-set redundency", value="gene_set_redundency_hist",
                                            plotOutput("gsc_boxplot2", height="340px", width="100%"),
                                            HTML("<div style='max-width: 600px'><b>Number of <font color='#9966ff'>gene-sets</font> per <font color='#009900'>gene</font>.</b>
                          Gives you a general overview of how many gene-sets each gene belongs to.
                          Do you have many genes contributing to a lot of overlap between gene-sets,
                          or do you have very uniquely defined gene-sets where genes only appear in one
                          or a few gene-sets?</div>")
                                   )
                            )
                     )
            ))),
            tabPanel(title="Gene-set table", value="tab_gsares",
                     div(id="loading_tab_gsares",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                     ),
                     hidden(div(id="tab_gsares",
                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                     h4(textOutput("gsaresTableText")),
                     HTML("</td><td>"),
                     uiOutput("gsaresTableText2"),
                     HTML("</td></tr></table>"),
                     DT::dataTableOutput('gsaresTable',height="100%"),
                     absolutePanel(id="controls", class="panel panel-default", fixed=FALSE,
                                   draggable=FALSE, top=110, left=20, right="auto", bottom="auto",
                                   width = "auto", height=70, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                   downloadButton("download_gs_table", "Download table", class="btn btn-primary btn-xs"),
                                   HTML("<span style='display:inline-block; width: 100px;'></span><b>Toggle columns:</b><span style='display:inline-block; width: 30px;'></span>"),
                                   HTML('Gene count<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_genes\', Math.random())" checked><span class="slider round"></span></label>'),
                                   HTML('Adjusted p-value<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_padj\', Math.random())" checked><span class="slider round"></span></label>'),
                                   HTML('p-value<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_p\', Math.random())"><span class="slider round"></span></label>'),
                                   HTML('Gene-set statistic<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_stats\', Math.random())"><span class="slider round"></span></label>')
                     )
                     ))
            ),
            tabPanel(title="Gene-set summary", value="tab_gssum",
                     div(id="loading_tab_gssum",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                     ),
                     hidden(div(id="tab_gssum",
                     fluidRow(column(12,
                                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                                     h4(textOutput("text_selected_gs")),
                                     HTML("</td><td>"),
                                     actionLink("link_tab_gsares", "Select another gene-set"),
                                     HTML("</td></tr></table>")
                     )),
                     fluidRow(column(5,
                                     fluidRow(column(12,
                                                     box(title="Gene-set stats", status="primary", width="100%",
                                                               tableOutput('gsTable2')
                                                     ),
                                                     box(title="Genes", status="primary", width="100%",
                                                         htmlOutput('gsNgenes'),
                                                         htmlOutput('link_view_gene'),
                                                         uiOutput("links_gsGenes")
                                                     ),
                                                     box(title="Gene-set neighbors", status="primary", width="100%",
                                                               HTML("Number of overlapping genes:"),
                                                               uiOutput('gsNeighborsSlider'),
                                                               actionLink("filter_gs2", "View gene-sets in table"),
                                                               HTML("(also makes gene-sets available in network plot)"),
                                                               uiOutput('gsNeighbors')
                                                     )
                                            )
                                     )
                            ),
                     column(7,
                            tabBox(width="100%",
                            #wellPanel(style="background-color: #ffffff; overflow-y:scroll; max-height: 1200px",
                                #      tabsetPanel(id="tabset_2", type="tabs",
                                                  tabPanel(title="Histogram", value="tab_histogram",
                                      HTML("<center><b>Histogram of all gene-level statistics,<br>"),
                                      HTML("and indiviual values for genes in selected gene-set.</b></center>"),
                                      plotOutput("gsHist", click=clickOpts(id="gsHist_click")),
                                      htmlOutput("geneStatType"), # just needed so that output.geneStatType is set... hacky
                                      conditionalPanel(condition="output.geneStatType == '<em></em>'", plotOutput("gsHist2", click=clickOpts(id="gsHist_click2")))
                                                  ),
                                      tabPanel(title="Boxplot", value="tab_boxplots",
                                      HTML("<center><b>Boxplots of gene-level statistics</b></center>"),
                                      plotOutput("gsBoxplots", click=clickOpts(id="gsBoxplots_click"))
                                      #plotOutput("gsRunningSum")
                                 #     )
                                      )
                            )
                     ))
                     ))
            ),
            tabPanel(title="Gene table", value="tab_genetable",
                     div(id="loading_tab_genetable",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                     ),
                     hidden(div(id="tab_genetable",
                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                     h4(textOutput("geneTableText")),
                     HTML("</td><td>"),
                     uiOutput("geneTableText2"),
                     HTML("</td></tr></table>"),
                     DT::dataTableOutput('geneTable', height="100%"),
                     absolutePanel(id="controls", class="panel panel-default", fixed=FALSE,
                                   draggable=FALSE, top=110, left=20, right="auto", bottom="auto",
                                   width = "auto", height=70, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                   downloadButton("download_gene_table", "Download table", class="btn btn-primary btn-xs"),
                                   HTML("<span style='display:inline-block; width: 100px;'></span>"),
                                   HTML('<b>Truncate long: </b><label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'truncate_gene_table\', Math.random())" checked><span class="slider round"></span></label>'),
                                   HTML("<span style='display:inline-block; width: 5px;'></span>"),
                                   HTML("(Hold pointer over a cell to see full content)")
                     )
                     ))
            ),
            tabPanel(title="Gene summary", value="tab_geneinfo",
                     div(id="loading_tab_geneinfo",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                     ),
                     hidden(div(id="tab_geneinfo",
                     fluidRow(column(12,
                                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                                     h4(textOutput("text_selected_gene")),
                                     HTML("</td><td>"),
                                     actionLink("link_tab_genetable", "Select another gene"),
                                     HTML("</td></tr></table>")
                     )),
                     fluidRow(column(6,
                                     box(title="Gene stats", status="primary", width="100%",
                                               htmlOutput("info_selected_gene")

                                     ),
                                     box(title="Gene-set info", status="primary", width="100%",
                                               htmlOutput("text_gene_gene_sets1"),
                                               actionLink("filter_gs", "View gene-sets in table"),
                                               HTML("(also makes gene-sets available in network plot)"),
                                               htmlOutput("text_gene_gene_sets2")
                                     )
                     ),
                     column(6,
                            box(title="Additional annotation", status="primary", width="100%",
                                      htmlOutput("anno_selected_gene")
                            )
                     ))
                     ))
            ),
            tabPanel(title="Network plot", value="tab_nwplot",
                     div(id="loading_tab_nwplot",
                         HTML("<br><br><center><h3>",c("Just hold on a sec","Please wait","Just a moment")[sample(1:3,1)],"</h3><h4>Loading your results...</h4></center>")
                     ),
                     hidden(div(id="tab_nwplot",
                     HTML("<br>"),
                     fluidRow(column(8,
                                     box(title=NULL, status="primary", height="85vh", width="100%", solidHeader =TRUE,
                                                visNetworkOutput("network", width = "100%", height = "80vh")
                                     ),
                                     absolutePanel(id="nw_geneset_table", class="panel panel-default", fixed=FALSE,
                                                   draggable=FALSE, top=10, left=30, right="auto", bottom="auto",
                                                   width=200, height=80, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                                   actionLink("filter_gs3", "View gene-sets in table")
                                     ),
                                     absolutePanel(id="nw_color_legend", class="panel panel-default", fixed=FALSE,
                                                   draggable=FALSE, top=10, left="auto", right=25, bottom="auto",
                                                   width=200, height=80, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                                   plotOutput("network_color_legend", width="200px", height="70px")
                                     )
                              ),
                              column(4,
                                     #wellPanel(style="background-color: #ffffff",
                                     box(title="Gene-set p-value to display", status="primary", collapsible=TRUE, width="100%",
                                               #fluidRow(column(12,HTML("<h4>Gene-set p-value to display</h4>"))),
                                               fluidRow(column(6,
                                                        selectInput("network_class", HTML("<b>P-value class</b>"),
                                                                     choices = selectable_classes_network,
                                                                     selected = selected_class_network)
                                                        ),
                                                        column(6,
                                                        selectInput("network_adjusted", HTML("<b>P-value adjustment</b>"),
                                                                     choices = list("Adjusted p-value" = TRUE,
                                                                                    "Non-adjusted p-value" = FALSE),
                                                                     selected = TRUE)
                                                        )
                                               )
                                     ),
                                     #wellPanel(style="background-color: #ffffff",
                                     box(title="Gene-set selection", status="primary", collapsible=TRUE, width="100%",
                                               #fluidRow(column(12,HTML("<h4>Gene-set selection</h4>"))),
                                               fluidRow(column(6,
                                                        selectInput("geneset_selection","Select gene-set by",
                                                                    choices=list("Significance cutoff"="significance",
                                                                                 "Predefined list"="list"),
                                                                    selected="significance")
                                               ),
                                               column(6,
                                                      conditionalPanel('input.geneset_selection=="significance"',
                                                      numericInput("network_significance", "Significance cutoff", cutoff_network,
                                                                   min=0, max=1, step=0.001)),
                                                      conditionalPanel('input.geneset_selection=="list"',
                                                                       selectInput("network_genesetlist", label="Select gene-sets from list",
                                                                  choices=setNames(as.list(names(genesets)),names(genesets))))
                                               )
                                               ),
                                               fluidRow(column(12,HTML("<b>Maximum allowed nodes</b>"),
                                                        numericInput("maxAllowedNodes", NULL, 50,
                                                                     min=0, max=Inf, step=1, width="100px")
                                                        )
                                               )
                                     ),
                                     #wellPanel(style="background-color: #ffffff",
                                     box(title="Network layout", status="primary", collapsible=TRUE, width="100%",
                                               #fluidRow(column(12,HTML("<h4>Network layout</h4>"))),
                                               fluidRow(column(8,
                                                        selectInput("layout", label=NULL,
                                                                    choices=list("Default (visNetwork)"="visNetwork",
                                                                                 "Nicely (igraph)"="layout_nicely",
                                                                                 "Fruchterman-Reingold (igraph)"="layout_with_fr",
                                                                                 "Kamada-Kawai (igraph)"="layout_with_kk",
                                                                                 "Sugiyama (igraph)"="layout_with_sugiyama",
                                                                                 "Star (igraph)"="layout_as_star",
                                                                                 "Circle (igraph)"="layout_in_circle",
                                                                                 "Grid (igraph)"="layout_on_grid"
                                                                                 ),
                                                                    selected="visNetwork")
                                                        ),
                                                        column(4,
                                                        conditionalPanel("input.layout == 'visNetwork' | input.layout == 'layout_nicely' | input.layout == 'layout_with_fr'",
                                                                                actionButton("layout_seed",label="Generate new layout", class="btn btn-primary btn-xs")
                                                        )
                                                        )
                                               ),
                                               fluidRow(column(12,
                                                        HTML('<b>Physics simulation: </b><label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'physics\', Math.random())" checked><span class="slider round"></span></label>')
                                               )
                                               )
                                     ),
                                     #wellPanel(style="background-color: #ffffff",
                                     box(title="Node properties", status="primary", collapsible=TRUE, width="100%", collapsed=TRUE,
                                               #fluidRow(column(12,HTML("<h4>Node properties</h4>"))),
                                               fluidRow(column(12,
                                                        sliderInput("nodeSize", label="Size range", min=1, max=100, value=c(10,40), ticks=FALSE)
                                                        )
                                               ),
                                               fluidRow(column(3,
                                                        selectInput("labelSize", label="Label size",
                                                                    choices=list("Off"=0,
                                                                                 "8"=8,
                                                                                 "10"=10,
                                                                                 "14"=14,
                                                                                 "18"=18,
                                                                                 "22"=22,
                                                                                 "30"=30,
                                                                                 "40"=40),
                                                                    selected=22)
                                                        ),
                                                        column(9,HTML("<br>"),
                                                        checkboxInput("truncate_labels", label="Truncate long labels", FALSE)
                                                        )
                                               )
                                     ),
                                     #wellPanel(style="background-color: #ffffff",
                                     box(title="Edge properties", status="primary", collapsible=TRUE, width="100%", collapsed=TRUE,
                                               #fluidRow(column(12,HTML("<h4>Edge properties</h4>"))),
                                               fluidRow(column(12,
                                                        sliderInput("edgeWidth", label="Width range", min=1, max=50, value=c(1,15), ticks=FALSE)
                                                        )
                                               ),
                                               fluidRow(column(12,HTML("<b>Required node overlap to draw edge</b>"))),
                                               fluidRow(column(3,
                                                        numericInput("edge_overlap", label=NULL, value=30, min=0, max=Inf)
                                                        ),
                                                        column(9,
                                                        selectInput("genes_or_percent", label=NULL,
                                                                    choices=list("genes"="genes",
                                                                                 "% of smallest gene-set in pair"="percent"),
                                                                    selected="percent")
                                                        )
                                               )
                                     )
                              )
                     )

                     ))
            ),
            tabPanel(title="Heatmap", value="tab_heatmap",
                     HTML("<br><br>This feature is currently unavailable but will be added soon! For now, use the <tt>GSAheatmap</tt> funtion in R.")
            )#,
            #tabPanel(title="Help", value="tab_help",
            #         HTML("<br><br>This feature is currently unavailable but will be added soon! For now, use the <tt>GSAheatmap</tt> funtion in R.")
            #)
          #)
        #)
      #),
    ))),

# ===================================================================================
    # server
    server <- function(input, output, session) {

      session$onSessionEnded(stopApp)

      # Output values used as parameters:
      output$geneStatType <- renderText(ifelse(gsares$geneStatType=="p-signed",'<em></em>','')) # super hacky...

      # Initialize parameters / static objects:
      gs_per_gene_list <- unstack(stack(gsares$gsc)[,2:1])
      gsares_directions <- switch(as.numeric(length(gsares$directions)==1)+1, gsares$directions, gsares$geneLevelStats*NA)
      colnames(gsares_directions) <- "directions"

      genetable_all <- stack(lapply(gs_per_gene_list,length))[,2:1]
      genetable_all <- merge(cbind(gsares$geneLevelStats, gsares_directions), genetable_all, by.x=0, by.y=1)
      rownames(genetable_all) <- genetable_all[,1]
      colnames(genetable_all) <- c("Gene ID","Gene-level statistic","Sign (FC direction)","In gene-sets")

      if(!is.null(geneAnnot)) {
        genetable_all <- merge(genetable_all, geneAnnot, by=1, all.x=TRUE)
        colnames(genetable_all) <- c("Gene ID","Gene-level statistic","Sign (FC direction)","In gene-sets", colnames(geneAnnot)[-1])
        rownames(genetable_all) <- genetable_all[,1]
      }
      if(all(is.na(genetable_all[,"Sign (FC direction)"]))) genetable_all <- genetable_all[,!colnames(genetable_all)%in%"Sign (FC direction)"]

      gsatab_colnames <- c("Name", "Genes (tot)", "Stat (dist.dir)", "Stat (dist.dir.up)", "p (dist.dir.up)", "p adj (dist.dir.up)",
                           "Stat (dist.dir.dn)", "p (dist.dir.dn)", "p adj (dist.dir.dn)", "Stat (non-dir.)", "p (non-dir.)", "p adj (non-dir.)",
                           "Genes (up)", "Stat (mix.dir.up)", "p (mix.dir.up)", "p adj (mix.dir.up)", "Genes (down)",
                           "Stat (mix.dir.dn)", "p (mix.dir.dn)", "p adj (mix.dir.dn)")
      gsatab_colnames <- gsatab_colnames[c(1,2,13,17,9,20,12,16,6,8,19,11,15,5,10,3,7,4,18,14)]
      gsatab_colnames <- gsatab_colnames[gsatab_colnames%in%colnames(GSAsummaryTable(gsares))]
      if(all(gsares$nGenesUp==0) & all(gsares$nGenesDn==0)) gsatab_colnames <- gsatab_colnames[!gsatab_colnames%in%c("Genes (up)","Genes (down)")]

      # Initialize reactive values:
      rval <- reactiveValues(sel_gs=names(gsares$gsc)[1],
                             sel_gene=gsares$gsc[[1]][1],
                             gs_table=GSAsummaryTable(gsares),
                             sel_gene_for_filtering=NULL,
                             sel_geneset_for_filtering=NULL,
                             red_gene="",
                             gene_table=genetable_all,
                             log_gs_plots=FALSE,
                             visible_columns=gsatab_colnames[!gsatab_colnames%in%c(grep("p \\(",gsatab_colnames, value=TRUE),grep("Stat",gsatab_colnames, value=TRUE))],
                             physics=TRUE,
                             colorLegendInfo=NULL,
                             nwGeneSets=NULL,
                             genesets=genesets,
                             maxNcharGenetable=TRUE)

      # -----------------------------------------------------------------------------
      #output$loggo <- renderImage({
      #  filename <- "loggo.png"
      #  list(src = filename,
      #       alt = "PIANO")
      #}, deleteFile = FALSE)


      # Run info --------------------------------------------------------------------
      output$text_runinfo1 <- renderUI({
        tmp <- gsares$info
        info <- c(
          paste("<b>Input:</b>",tmp$removedGSnoGenes + tmp$removedGSsizeLimit + tmp$nGeneSets,"gene-sets,",tmp$nGenesStatistics,"genes<br>"),
          paste("- Removed",tmp$removedGSsizeLimit,"gene-sets not matching the size limits<br>"),
          paste("- Removed",tmp$removedGenesGSC,"genes from GSC due to lack of matching gene-level statistics<br>"),
          paste("- Removed",tmp$removedGSnoGenes,"gene-sets containing no genes after gene removal<br>"),
          paste("- Out of the",tmp$nGenesStatistics,"genes,",tmp$nGenesGSC,"are associated with gene-sets<br>"),
          paste("<b>Final gene/gene-set association:</b>",tmp$nGeneSets,"gene-sets,",tmp$nGenesGSC,"genes*")
        )
        HTML(info)
      })
      output$text_runinfo2 <- renderUI({
        tmp <- gsares$info
        info <- c(
          paste("<b>Total run time:</b>",round(gsares$runtime[3]/60,2),"min<br>"),
          paste("<b>GSA method:</b>",gsares$geneSetStat,"<br>"),
          paste("<b>Gene-set statistic name:</b>",gsares$gsStatName,"<br>"),
          paste("<b>Gene-set size limits:</b> min",paste(gsares$gsSizeLim,collapse=", max "),"genes<br>"),
          paste("<b>Significance calculation:</b>",gsares$signifMethod,"<br>"),
          paste("<b>Permutations:</b> ",gsares$nPerm,"*<br>",sep=""),
          paste("<b>P-value adjustment:</b>",gsares$adjMethod)
        )
        HTML(info)
      })
      output$text_runinfo3 <- renderUI({
          HTML("<em>*) In case premuted gene-level statistics are used for significance estimation, all input gene-level
          statistics are used. Permutations are also used by the reporter features method to calculate gene-set statistics,
          i.e. regardless of significance estimation method.</em>")
      })

      output$gsc_boxplot1 <- renderPlot({
        par(mfcol=c(2,1))
        layout(mat = matrix(c(1,2),2,1, byrow=TRUE), heights= c(3,2))
        tmp <- unlist(lapply(gsares$gsc,length))
        par(mar=c(0, 4.1, 1.1, 2.1))
        h <- hist(tmp, n=50, xlim=c(0,max(tmp)), xaxt="n", xlab=NULL, main=NULL, col="forestgreen")
        par(mar=c(5, 4.1, 0, 2.1))
        boxplot(tmp, horizontal=TRUE, ylim=c(0,max(tmp)), col="forestgreen", axes=FALSE, xlab="Number of genes")
        axis(1, las=2)
      })
      output$gsc_boxplot2 <- renderPlot({
        layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  heights = c(3,2))
        tmp <- unlist(lapply(gs_per_gene_list,length))
        par(mar=c(0, 4.1, 1.1, 2.1))
        h <- hist(tmp, n=50, xlim=c(0,max(tmp)), xaxt="n", xlab=NULL, main=NULL, col="mediumpurple2")
        par(mar=c(5, 4.1, 0, 2.1))
        boxplot(tmp, horizontal=TRUE, ylim=c(0,max(tmp)), col="mediumpurple2", axes=FALSE, xlab="Number of gene-sets")
        axis(1, las=2)
      })

      # GSAres table ----------------------------------------------------------------
      output$gsaresTableText <- renderText({
         if(nrow(rval$gs_table)<length(gsares$gsc)) {
           rval$geneset_filtering_text
         } else {
           "Showing all gene-sets"
         }
      })
      output$gsaresTableText2 <- renderUI({
        if(nrow(rval$gs_table)<length(gsares$gsc)) {
          actionLink("filter_all_gs","(Show all gene-sets)")
        }
      })
      observeEvent(input$filter_all_gs, {
        rval$sel_gene_for_filtering <- NULL
        rval$gs_table <- GSAsummaryTable(gsares)
      })

      # Toggle columns events:
      observeEvent(input$toggle_genes, {
        tmp <- grep("Genes",gsatab_colnames, value=TRUE)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_padj, {
        tmp <- grep("p adj",gsatab_colnames, value=TRUE)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_p, {
        tmp <- grep("p \\(",gsatab_colnames, value=TRUE)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_stats, {
        tmp <- grep("Stat",gsatab_colnames, value=TRUE)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })

      output$gsaresTable <- DT::renderDataTable({tmp <- rval$gs_table[gsatab_colnames[gsatab_colnames%in%rval$visible_columns]]},
                                                server=TRUE, escape=FALSE,
                                                selection=list(mode='single', target='row'),
                                                filter="none",
                                                rownames=FALSE,
                                                fillContainer=TRUE,
                                                extensions=c("Scroller"),
                                                options = list(initComplete=JS(
                                                    # this is to set color format of first row
                                                    "function(settings, json) {",
                                                    "$(this.api().table().header()).css({'background-color': '#E2F3FF', 'color': '#000'});",
                                                    "}"),
                                                    dom = 'frtip',
                                                    deferRender=TRUE,
                                                    scrollY="calc(100vh - 250px)",
                                                    scrollX=TRUE,
                                                    scroller=TRUE
                                                ))

      observeEvent(input$gsaresTable_rows_selected, {
        rval$sel_gs <- rval$gs_table$Name[input$gsaresTable_rows_selected]
        rval$red_gene <- ""
        updateNavbarPage(session,"navbarpage",selected="tab_gssum")
        selectRows(dataTableProxy("gsaresTable", session), list())
      })
      
      output$download_gs_table <- downloadHandler(
        filename = function() {
          paste("gene-set_table", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(rval$gs_table[gsatab_colnames[gsatab_colnames%in%rval$visible_columns]], file, row.names=FALSE)
        }
      )
      

      # Gene-set summary ------------------------------------------------------------
      output$text_selected_gs <- renderText({
        rval$sel_gs
      })

      observeEvent(input$link_tab_gsares, {
        updateNavbarPage(session,"navbarpage",selected="tab_gsares")
      })

      # Tables with info:
      output$gsTable2 <- renderTable({
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        tmp <- tmp$stats
        rownames(tmp) <- tmp[,1]
        tmp <- rbind(
        format(rev(c(tmp["p (dist.dir.dn)",2],tmp["p (mix.dir.dn)",2],tmp["p (non-dir)",2],tmp["p (mix.dir.up)",2],tmp["p (dist.dir.up)",2])), digits=3, scientific=TRUE),
        format(rev(c(tmp["p adj (dist.dir.dn)",2],tmp["p adj (mix.dir.dn)",2],tmp["p adj (non-dir)",2],tmp["p adj (mix.dir.up)",2],tmp["p adj (dist.dir.up)",2])), digits=3, scientific=TRUE)
        )
        rownames(tmp) <- c("p-value","Adjusted p-value")
        colnames(tmp) <- rev(c("Distinct directional (down)","Mixed directional (down)","Non-directional","Mixed directional (up)", "Distinct directional (up)"))
        t(tmp)
      }, rownames=TRUE,colnames=TRUE)

      output$gsNgenes <- renderUI({
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        tmp_tot <- tmp$stats[tmp$stats$Name=="Genes (tot)",2]
        tmp_up <- tmp$stats[tmp$stats$Name=="Genes (up)",2]
        tmp_dn <- tmp$stats[tmp$stats$Name=="Genes (dn)",2]
        if(tmp_up>0) {
          tmp_up <- paste("up:",tmp_up)
        } else {
          tmp_up <- NULL
        }
        if(tmp_dn>0) {
          tmp_dn <- paste("down:",tmp_dn)
        } else {
          tmp_dn <- NULL
        }
        tmp <- ""
        if(!is.null(tmp_up) | !is.null(tmp_dn)) tmp <- paste("(",paste(tmp_up, tmp_dn, collapse=", )"),")",sep="")
        tmp <- paste("<b>The gene-set contains ",tmp_tot,
                     " genes </b>", tmp, "<br>",
                     sep=""
        )
        HTML(tmp)
      })

      # Genes in gene-set info:
      output$link_view_gene <- renderUI({
        if(rval$red_gene!="") {
          list(
            HTML(paste("Selected:",rval$red_gene)),
            actionLink("view_gene","(view details)"),
            HTML("<br><br>"),
            actionLink("view_genes_table",paste("View genes in table"))
          )
        } else if(rval$red_gene=="") {
          list(
            HTML("<em><font color='#ff0000'>You can click genes below, or in the plots to the right...</font></em><br><br>"),
            actionLink("view_genes_table",paste("View genes in table"))
          )
        }
      })
      observeEvent(input$view_gene, {
        rval$sel_gene <- rval$red_gene
        updateNavbarPage(session,"navbarpage",selected="tab_geneinfo")
      })
      observeEvent(input$view_genes_table, {
        rval$gene_table <- genetable_all[genetable_all$`Gene ID` %in% names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats),]
        rval$sel_geneset_for_filtering <- rval$sel_gs
        updateNavbarPage(session,"navbarpage",selected="tab_genetable")
      })

      # Genes in gene-set list:
      output$links_gsGenes <- renderUI({
        tmp <- geneSetSummary(gsares,rval$sel_gs)$geneLevelStats
        tmp <- sort(tmp)
        tmp <- names(tmp)

        tmp1 <- ifelse(tmp == rval$red_gene,"<b><font color='#ff0000'>","<font color='#000000'>")
        tmp2 <- ifelse(tmp == rval$red_gene,"</font></b>","</font>")

        HTML(paste(paste("<a href='#'"," onclick='Shiny.onInputChange(",'"links_gsGenes_click", "',tmp,'"',");'",">",tmp1,tmp,tmp2,"</a>",sep=""), collapse=", "))
      })
      observeEvent(input$links_gsGenes_click, {
        rval$red_gene <- input$links_gsGenes_click
      })

      # Neighbor gene-sets:
      output$gsNeighborsSlider <- renderUI({
        tmp <- length(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        sliderInput("nNeighborGenesets", NULL, 1, tmp, min(c(5,tmp)), round=TRUE, step=1, ticks=FALSE, width="300px")
      })
      output$gsNeighbors <- renderUI({
        tmp <- names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        tmp <- setdiff(names(gsares$gsc)[unlist(lapply(gsares$gsc, function(x) sum(tmp%in%x)>=input$nNeighborGenesets))], rval$sel_gs)
        HTML(paste("<b>",rval$sel_gs,"</b>,",sep=""),paste(paste("<a href='#'"," onclick='Shiny.onInputChange(",'"links_genesets_click", "',tmp,'"',");'","><font color='#000000'>",tmp,"</font></a>",sep=""),collapse=", "))
      })
      observeEvent(input$filter_gs2, { # action link in ui
        tmp1 <- names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        tmp1 <- setdiff(names(gsares$gsc)[unlist(lapply(gsares$gsc, function(x) sum(tmp1%in%x)>=input$nNeighborGenesets))], rval$sel_gs)
        tmp2 <- GSAsummaryTable(gsares)
        rval$gs_table <- tmp2[tmp2$Name %in% c(tmp1,rval$sel_gs),]
        rval$geneset_filtering_text <- paste("Gene-sets neighboring",rval$sel_gs)
        rval$genesets <- setNames(append(rval$genesets,list(rval$gs_table$Name)),c(names(rval$genesets),rval$geneset_filtering_text))
        if(any(names(rval$genesets) == "No available gene-sets")) {
          rval$genesets <- rval$genesets[rval$genesets != "No available gene-sets"]
        }
        updateSelectInput(session, "network_genesetlist",
                          choices=setNames(names(rval$genesets),names(rval$genesets)))
        updateNavbarPage(session,"navbarpage",selected="tab_gsares")
      })

      # Plots:

      # Histogram:
      output$gsHist <- renderPlot({

        par(mar=c(5,10,1,5))
        tmp2 <- gsares$geneLevelStats
        if(gsares$geneStatType=="p-signed") tmp2 <- tmp2[gsares$directions>0]
        hist(tmp2,100,col="white", xlab=paste("Gene-level statistics",ifelse(gsares$geneStatType=="p-signed"," (direction: up)",""),sep=""), main="")
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        if(gsares$geneStatType=="p-signed") {
          tmp <- tmp$geneLevelStats[tmp$directions>0]
        } else {
          tmp <- tmp$geneLevelStats
        }
        abline(v=tmp, lty=2)
        hist(tmp2,100,col=ifelse(gsares$geneStatType=="p-signed","indianred2","grey"), add=TRUE)

        # Highlight selected gene:
        observeEvent(input$gsHist_click, {
          rval$red_gene <- names(tmp)[which.min(abs(tmp - input$gsHist_click$x))]
        })
        abline(v=tmp[rval$red_gene], col="red", lwd=2)
      })


      output$gsHist2 <- renderPlot({
        par(mar=c(5,10,1,5))
        tmp2 <- gsares$geneLevelStats[gsares$directions<=0]
        hist(tmp2,100,col="white", xlab="Gene-level statistics (direction: down)", main="")
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        tmp <- tmp$geneLevelStats[tmp$directions<=0]
        abline(v=tmp, lty=2)
        hist(tmp2,100,col="skyblue2", add=TRUE)

        # Highlight selected gene:
        observeEvent(input$gsHist_click2, {
          rval$red_gene <- names(tmp)[which.min(abs(tmp - input$gsHist_click2$x))]
        })
        abline(v=tmp[rval$red_gene], col="dodgerblue3", lwd=2)
      })

      # Boxplots:
      output$gsBoxplots <- renderPlot({

        par(mar=c(5,10,1,5))
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        gls_gs <- tmp$geneLevelStats
        gls_gs_up <- gls_gs[tmp$directions>0]
        gls_gs_dn <- gls_gs[tmp$directions<=0]
        gls_all <- gsares$geneLevelStats
        gls_all_up <- gls_all[gsares$directions>0]
        gls_all_dn <- gls_all[gsares$directions<=0]

        if(gsares$geneStatType=="p-signed") {
          bp <- boxplot(list(gls_all_dn, gls_all_up, gls_gs_dn, gls_gs_up),
                        horizontal=TRUE,
                        names=c(paste("All genes down (",length(gls_all_dn),")",sep=""),
                                paste("All genes up (",length(gls_all_up),")",sep=""),
                                "Gene-set genes (down)","Gene-set genes (up)"),
                        las=1,
                        col=c("skyblue2","indianred2","white","white"),
                        border=c("black","black","dodgerblue3","red"),
                        main="", xlab="Gene-level statistics",
                        outline=TRUE, outcol="white")

          # Add outliers:
          tmp <- cbind(bp$out,bp$group)
          tmp <- tmp[tmp[,2] %in% 1:2,]
          points(tmp)

          # Add jitter:
          set.seed(1)
          y_dn <- jitter(rep(3,length(gls_gs_dn)), amount=0.3)
          names(y_dn) <- names(gls_gs_dn)
          points(x=gls_gs_dn, y=y_dn)

          y_up <- jitter(rep(4,length(gls_gs_up)), amount=0.3)
          names(y_up) <- names(gls_gs_up)
          points(x=gls_gs_up, y=y_up)

        } else {
          bp <- boxplot(list(gls_all, gls_gs),
                        horizontal=TRUE, names=c(paste("All genes (",length(gls_all),")",sep=""),"Genes in gene-set"), las=1,
                        col=c("grey","white"),
                        main="", xlab="Gene-level statistics",
                        outline=TRUE, outcol="white")
          # Add outliers:
          tmp <- cbind(bp$out,bp$group)
          tmp <- tmp[tmp[,2] == 1,]
          points(tmp)

          # Add jitter:
          set.seed(1)
          y <- jitter(rep(2,length(gls_gs)), amount=0.3)
          names(y) <- names(gls_gs)
          points(x=gls_gs, y=y)
        }

        # Highlight selected point:
        if(gsares$geneStatType=="p-signed") {
          x <- c(gls_gs_up, gls_gs_dn)
          y <- c(y_up, y_dn)
        } else {
          x <- gls_gs
        }
        observeEvent(input$gsBoxplots_click, {
          rval$red_gene <- names(y)[which.min(as.matrix(dist(rbind(input$gsBoxplots_click[c("x","y")],cbind(x,y))))[-1,1])]
        })
        points(x=x[rval$red_gene], y=y[rval$red_gene], pch=16, cex=2, col=ifelse(rval$red_gene%in%names(gls_gs_up),"red","dodgerblue3"))
      })

      # Running sum plot:
      #output$gsRunningSum <- renderPlot({
      #  par(mar=c(3,10,4.1,5))
      #  plotRunningSum(gsares, rval$sel_gs, gseaParam=1)
      #})


      # Gene table ----------------------------------------------------------

      output$geneTableText <- renderText({
        if(nrow(rval$gene_table) < nrow(genetable_all)) {
          paste("Genes in ", rval$sel_geneset_for_filtering)
        } else {
          "Showing all genes in GSC"
        }
      })
      output$geneTableText2 <- renderUI({
        if(nrow(rval$gene_table) < nrow(genetable_all)) {
          actionLink("filter_all_genes","(Show all genes)")
        }
      })
      observeEvent(input$filter_all_genes, {
        rval$sel_geneset_for_filtering <- NULL
        rval$gene_table <- genetable_all
      })

      output$geneTable <- DT::renderDataTable({apply(cbind(rank(rval$gene_table$`Gene-level statistic`),rval$gene_table),2,function(x) {
                                                  if(rval$maxNcharGenetable) {
                                                    tmp <- 50
                                                  } else {
                                                    tmp <- 10000
                                                  }
                                                  paste("<div title='",ifelse(is.na(x),"",x),
                                                        "' class='DT_genetable_cell'><p>",
                                                        ifelse(is.na(x),"",substr(x,1,tmp)),
                                                        ifelse(is.na(x),"",ifelse(nchar(x)>tmp,"...","")),
                                                        "</div>",sep="")
                                                }
                                                )
                                               },
        server=TRUE, escape=FALSE,
        selection=list(mode='single', target='row'),
        filter="none",
        rownames=FALSE,
        extensions=c("Scroller"),#,"FixedColumns"),
        fillContainer=TRUE,
        options = list(initComplete=JS(
        # this is to set color format of first row
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#E2F3FF', 'color': '#000'});",
        "}"),
        dom = 'frtip',
        deferRender=TRUE,
        scrollY="calc(100vh - 250px)",
        scrollX=TRUE,
        scroller=TRUE,
        columnDefs=list(list(orderData=0, targets=2),
                        list(visible=FALSE, targets=0))
        #fixedColumns=list(leftColumns=1))
        )
      )

      observeEvent(input$geneTable_rows_selected, {
        print(class(rval$gene_table[,"Gene-level statistic"]))
        rval$sel_gene <- rval$gene_table$`Gene ID`[input$geneTable_rows_selected]
        rval$red_gene <- ""
        updateNavbarPage(session,"navbarpage",selected="tab_geneinfo")
        selectRows(dataTableProxy("geneTable", session), list())
      })

      observeEvent(input$truncate_gene_table, {
        rval$maxNcharGenetable <- ifelse(rval$maxNcharGenetable,FALSE,TRUE)
      })
      
      output$download_gene_table <- downloadHandler(
        filename = function() {
          paste("gene_table", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(rval$gene_table, file, row.names=FALSE)
        }
      )


      # Gene summary --------------------------------------------------------
      observeEvent(input$link_tab_genetable, { # action link in ui
        updateNavbarPage(session,"navbarpage",selected="tab_genetable")
      })

      output$text_selected_gene <- renderText({
        rval$sel_gene
      })

      output$info_selected_gene <- renderUI({
        HTML(paste("<b>Gene-level statistic:</b>",gsares$geneLevelStats[rval$sel_gene,],"<br>",
                   "<b>Direction:</b>",ifelse(is.na(gsares_directions[rval$sel_gene,]),"NA",ifelse(gsares_directions[rval$sel_gene,]>0,"Up","Down")))
        )
      })

      output$anno_selected_gene <- renderUI({
        HTML(paste(ifelse(ncol(genetable_all)>4,paste("<b>",colnames(genetable_all)[-c(1,2,3,4)], ":</b>",t(genetable_all[rval$sel_gene,-c(1,2,3,4)]), collapse="<br>"),""))
        )
      })

      output$text_gene_gene_sets1 <- renderUI({HTML(paste("<b>The gene is present in",genetable_all[rval$sel_gene,"In gene-sets"],"gene-sets:</b>"))})

      output$text_gene_gene_sets2 <- renderUI({
        tmp <- gs_per_gene_list[[rval$sel_gene]]
        HTML(paste(paste("<a href='#'"," onclick='Shiny.onInputChange(",'"links_genesets_click", "',tmp,'"',");'","><font color='#000000'>",tmp,"</font></a>",sep=""),collapse=", "))
      })
      observeEvent(input$filter_gs, { # action link in ui
        tmp <- GSAsummaryTable(gsares)
        rval$sel_gene_for_filtering <- rval$sel_gene
        rval$geneset_filtering_text <- paste("Gene-sets containing", rval$sel_gene_for_filtering)
        rval$gs_table <- tmp[tmp$Name %in% gs_per_gene_list[[rval$sel_gene]],]

        rval$genesets <- setNames(append(rval$genesets,list(rval$gs_table$Name)),c(names(rval$genesets),rval$geneset_filtering_text))
        if(any(names(rval$genesets) == "No available gene-sets")) {
          rval$genesets <- rval$genesets[rval$genesets != "No available gene-sets"]
        }
        updateSelectInput(session, "network_genesetlist",
                          choices=setNames(names(rval$genesets),names(rval$genesets)))

        updateNavbarPage(session,"navbarpage",selected="tab_gsares")

      })
      observeEvent(input$links_genesets_click, {
        rval$sel_gs <- input$links_genesets_click
        rval$red_gene <- ""
        selectRows(dataTableProxy("gsaresTable", session), list())
        updateNavbarPage(session,"navbarpage",selected="tab_gssum")
      })

      # Network plot --------------------------------------------------------------------

      output$network <- renderVisNetwork({

        if(input$geneset_selection=="significance") {
          network_significance <- input$network_significance
          network_genesetlist <- NULL
        } else if(input$geneset_selection=="list") {
          network_significance <- 1
          network_genesetlist <- rval$genesets[[input$network_genesetlist]]
        }

        nwPlot <- try(suppressWarnings(networkPlot2(gsares,
                               class=unlist(strsplit(input$network_class,"_"))[1],
                               direction=unlist(strsplit(input$network_class,"_"))[2],
                               adjusted=as.logical(input$network_adjusted),
                               significance=network_significance,
                               geneSets=network_genesetlist,
                               nodeSize=input$nodeSize,
                               labelSize=input$labelSize,
                               ncharLabel=ifelse(input$truncate_labels,10,Inf),
                               lay=input$layout,
                               physics=rval$physics,
                               seed=input$layout_seed,
                               edgeWidth=input$edgeWidth,
                               overlap=switch(input$genes_or_percent, genes=input$edge_overlap, percent=input$edge_overlap/100),
                               maxAllowedNodes=input$maxAllowedNodes,
                               shiny=TRUE,
                               main=NULL,
                               submain=NULL)), silent=TRUE)
        if(is(nwPlot[1], "try-error")) {
          nw_notification_text <- gsub(".*: (.*)","\\1",nwPlot[1])
          if(grepl("less than two gene sets were selected, can not plot",nwPlot[1])) {
            nw_notification_text <- "For the given parameters, less than two gene sets were selected. Try adjusting the significance cutoff."
          } else if(grepl("the selected parameters results in a network with more than",nwPlot[1])) {
            nw_notification_text <- paste("The selected parameters results in a network with more than",input$maxAllowedNodes,
                                          "nodes (gene-sets). Drawing large networks requires more memory, if you want to continue, increase the value of 'Maximum allowed nodes'",
                                          "or adjust the parameters to generate a smaller network, e.g. by decreasing the 'Significance cutoff'.",
                                          "Tip: Turn off the 'Physics simulation' for faster drawing of large networks.")
          } else if(grepl("argument geneSets not matching gene-set names in argument gsaRes",nwPlot[1])) {
            nw_notification_text <- "No gene-set lists available. Or all available lists contain incorrect gene-set names."
          }

          showNotification(nw_notification_text,
                           id="nw_notification",
                           duration=NULL,
                           closeButton=TRUE,
                           type="error")
          rval$colorLegendInfo <- NULL # to fail (= hide) the colorlegend plot
          cat() # needed to avoid visNetwork error print in console
        } else {
          removeNotification("nw_notification")
          rval$colorLegendInfo <- nwPlot$colorLegendInfo
          rval$nwGeneSets <- nwPlot$x$nodes$geneSetNames
          suppressWarnings(nwPlot %>% visExport(type="png", name="gene-set-network", 
                                                style="color: #ffffff; background-color: #3F8CBA; border-color: #2780e3; font-size: 13px;"))
        }
      })

      observeEvent(input$navbarpage, {
        if(input$navbarpage != "tab_nwplot")
          removeNotification("nw_notification")
      })

      observeEvent(input$physics, {
        rval$physics <- ifelse(rval$physics,FALSE,TRUE)
      })

      output$network_color_legend <- renderPlot({
        plotColorLegend <- function() {
          par(mar=c(3.5,0,0,0), oma=c(0,0,0,0))
          plot(1,1,col="white",cex=0, bty="n", ylab="", xlab="", yaxt="n", ylim=c(0,1), xlim=rval$colorLegendInfo$range)
          mtext("-log10(p)", side=1, line=2)
          xleft <- seq(from=rval$colorLegendInfo$range[1], to=rval$colorLegendInfo$range[2], length.out=length(rval$colorLegendInfo$colors))
          xright <- xleft + abs(xleft[1]-xleft[2])
          rect(xleft,0,xright,1, col=rval$colorLegendInfo$colors, border=NA)
        }
        tmp <- try(plotColorLegend(), silent=TRUE)
        if(is(tmp[1], "try-error")) plot(1,1,col="white",cex=0, bty="n", ylab="", xlab="", axes=FALSE)
      })

      observeEvent(input$filter_gs3, { # action link in ui
        tmp <- GSAsummaryTable(gsares)
        rval$geneset_filtering_text <- paste("Gene-sets in the network plot")
        rval$gs_table <- tmp[tmp$Name %in% rval$nwGeneSets,]
        updateNavbarPage(session,"navbarpage",selected="tab_gsares")
      })


      # Heatmap -------------------------------------------------------------------------
      #output$heatmap <- renderPlot({
      #  GSAheatmap(gsares, ncharLabel=400, cex=2)
      #})

      hide("loading_tab_gsainfo")
      showElement("tab_gsainfo")
      hide("loading_tab_gsares")
      showElement("tab_gsares")
      hide("loading_tab_gssum")
      showElement("tab_gssum")
      hide("loading_tab_genetable")
      showElement("tab_genetable")
      hide("loading_tab_geneinfo")
      showElement("tab_geneinfo")
      hide("loading_tab_nwplot")
      showElement("tab_nwplot")


    }
  )

  # launch app
  if(!browser) browser <- getOption("shiny.launch.browser", interactive())
  runApp(app, launch.browser=browser)

}
