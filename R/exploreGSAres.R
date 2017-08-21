exploreGSAres <- function(gsares, browser=T, geneAnnot=NULL) {
  
  # App object with ui and server
  app <- shinyApp(
    
# ===================================================================================    
    # ui
    ui <- fluidPage(
      includeCSS(system.file("shiny", "css", "cosmo.css", package="piano")),
      includeCSS(system.file("shiny", "css", "toggle_switch.css", package="piano")),
      includeCSS(system.file("shiny", "css", "style.css", package="piano")),
      
      fluidRow(
        column(8,
               h2(HTML("Explore gene-set analysis results"))
               ),
        column(4,
               div(imageOutput(system.file("shiny", "loggo.png", package="piano"), height="70px"), style="text-align: right;")
               #
        )
      ),
      fluidRow(
        column(12,
          tabsetPanel(id="tabset1",
            type = "tabs", 
            tabPanel(title="Run info", value="tab_gsainfo",
                     HTML("<br>"),
                     column(5,
                            wellPanel(style="background-color: #ffffff;",
                                      htmlOutput("text_runinfo1")
                            ),
                            wellPanel(style="background-color: #ffffff;",
                                      htmlOutput("text_runinfo2")
                            ),
                            wellPanel(style="background-color: #E2F3FF;",
                                      htmlOutput("text_runinfo3")
                            )
                     ),
                     column(7,
                            wellPanel(style="background-color: #ffffff; max-width: 1200px; overflow-y:scroll; max-height: calc(100vh - 200px)",
                                      HTML("<h4>Summary plots</h4>"),
                                      plotOutput("gsc_boxplot1", height="340px", width="100%"),
                                      HTML("<div style='max-width: 600px'><b>Fig 1. Number of <font color='#009900'>genes</font> per <font color='#9966ff'>gene-set</font>.</b>
                          Gives you an idea about sizes of the gene-sets used in the analysis.
                          Are you looking at small specific gene-sets with few genes, or large
                          general gene-sets? Note that these stats are <em>after</em> filtering the
                          input GSC according to available gene-level statistics and the gsSizeLim parameter.
                          <br><br></div>"),
                                      plotOutput("gsc_boxplot2", height="340px", width="100%"),
                                      HTML("<div style='max-width: 600px'><b>Fig 2. Number of <font color='#9966ff'>gene-sets</font> per <font color='#009900'>gene</font>.</b> 
                          Gives you a general overview of how many gene-sets each gene belongs to.
                          Do you have many genes contributing to a lot of overlap between gene-sets, 
                          or do you have very uniquely defined gene-sets where genes only appear in one
                          or a few gene-sets?</div>")
                            )
                     )
            ),
            tabPanel(title="Gene-set table", value="tab_gsares",
                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                     h4(textOutput("gsaresTableText")),
                     HTML("</td><td>"),
                     uiOutput("gsaresTableText2"),
                     HTML("</td></tr></table>"),
                     DT::dataTableOutput('gsaresTable',height="100%"),
                     absolutePanel(id="controls", class="panel panel-default", fixed=TRUE,
                                   draggable=F, top=150, left=20, right="auto", bottom="auto",
                                   width = "auto", height=70, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                   actionButton("dowload", "Download table", class="btn btn-primary btn-xs"),
                                   HTML("<span style='display:inline-block; width: 100px;'></span><b>Toggle columns:</b><span style='display:inline-block; width: 30px;'></span>"),
                                   HTML('Gene count<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_genes\', Math.random())" checked><span class="slider round"></span></label>'),
                                   HTML('Adjusted p-value<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_padj\', Math.random())" checked><span class="slider round"></span></label>'),
                                   HTML('p-value<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_p\', Math.random())"><span class="slider round"></span></label>'),
                                   HTML('Gene-set statistic<label class="switch"><input type="checkbox" onclick="Shiny.onInputChange(\'toggle_stats\', Math.random())"><span class="slider round"></span></label>')
                     )
            ),
            tabPanel(title="Gene-set summary", value="tab_gssum",
                     fluidRow(column(12,
                                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                                     h4(textOutput("text_selected_gs")),
                                     HTML("</td><td>"),
                                     actionLink("link_tab_gsares", "Select another gene-set"),
                                     HTML("</td></tr></table>")
                     )),
                     fluidRow(column(5,
                                     wellPanel(style="background-color: #ffffff;",
                                               tableOutput('gsTable2')
                                     ),
                                     wellPanel(style="background-color: #ffffff; overflow-y:scroll; max-height: 400px",
                                       htmlOutput('gsNgenes'),
                                       htmlOutput('link_view_gene'),
                                       uiOutput("links_gsGenes")
                                     ),
                                     wellPanel(style="background-color: #ffffff; overflow-y:scroll; max-height: 400px",
                                               HTML("<b>Gene-set neighbors</b>"),
                                               HTML("<table width='100%'><tr><td>Number of overlapping genes:</td><td align='left'>"),
                                               uiOutput('gsNeighborsSlider'),
                                               HTML("</td></tr></table>"),
                                               actionLink("filter_gs2", "View gene-sets in table"),
                                               uiOutput('gsNeighbors')
                                     )
                     ),
                     column(7,
                            wellPanel(style="background-color: #ffffff; overflow-y:scroll; max-height: 1200px",
                                      HTML("<center><b>Histogram of all gene-level statistics,<br>"),
                                      HTML("and indiviual values for genes in selected gene-set.</b></center>"),
                                      plotOutput("gsHist", click=clickOpts(id="gsHist_click")),
                                      htmlOutput("geneStatType"), # just needed so that output.geneStatType is set... hacky
                                      conditionalPanel(condition="output.geneStatType == '<em></em>'", plotOutput("gsHist2", click=clickOpts(id="gsHist_click2"))),
                                      HTML("<center><b>Boxplots of gene-level statistics</b></center>"),
                                      plotOutput("gsBoxplots", click=clickOpts(id="gsBoxplots_click"))
                                      #plotOutput("gsRunningSum")
                            )
                     ))
            ),
            tabPanel(title="Gene table", value="tab_genetable",
                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                     h4(textOutput("geneTableText")),
                     HTML("</td><td>"),
                     uiOutput("geneTableText2"),
                     HTML("</td></tr></table>"),
                     DT::dataTableOutput('geneTable', height="100%"),
                     absolutePanel(id="controls", class="panel panel-default", fixed=TRUE,
                                   draggable=F, top=160, left=20, right="auto", bottom="auto",
                                   width = "auto", height=70, style="border: 0px; box-shadow: 0px 0px; background-color: rgba(255, 255, 255, 0);",
                                   HTML("<table><tr><td>"),
                                   actionButton("dowload", "Download table", class="btn btn-primary btn-xs"),
                                   HTML("<span style='display:inline-block; width: 100px;'></td><td><b>Max characters per column:</b><span style='display:inline-block; width: 10px;'></td><td>"),
                                   uiOutput("ncharSlider"),
                                   HTML("</td><td width='10px'></td><td>(Hold pointer over a cell to see full content)</td></tr></table>")
                     )
            ),
            tabPanel(title="Gene summary", value="tab_geneinfo",
                     fluidRow(column(12,
                                     HTML("<table style='border-collapse:separate; border-spacing:10px;'><tr><td>"),
                                     h4(textOutput("text_selected_gene")),
                                     HTML("</td><td>"),
                                     actionLink("link_tab_genetable", "Select another gene"),
                                     HTML("</td></tr></table>")                
                     )),
                     wellPanel(style="background-color: #ffffff; max-width: 700px",
                     htmlOutput("info_selected_gene")
                     ),
                     wellPanel(style="background-color: #ffffff; max-width: 700px",
                     htmlOutput("text_gene_gene_sets1"),
                     actionLink("filter_gs", "View gene-sets in table"),
                     htmlOutput("text_gene_gene_sets2")
                     )
            ),
            tabPanel(title="Network plot", value="tab_nwplot",
                     #plotOutput("nwPlot", width="100%", height="800px")
                     HTML("<br><br>This feautere is currently unavailable but will be added soon! For now, use the <tt>networkPlot</tt> funtion in R.")
            ),
            tabPanel(title="Heatmap", value="tab_heatmap",
                     #plotOutput("heatmap", width="100%", height="800px")
                     HTML("<br><br>This feautere is currently unavailable but will be added soon! For now, use the <tt>GSAheatmap</tt> funtion in R.")
            )
          )
        )
      )
    ), 

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
        genetable_all <- merge(genetable_all, geneAnnot, by=1, all.x=T)
        colnames(genetable_all) <- c("Gene ID","Gene-level statistic","Sign (FC direction)","In gene-sets", colnames(geneAnnot)[-1])
        rownames(genetable_all) <- genetable_all[,1]
      }
      
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
                             visible_columns=gsatab_colnames[!gsatab_colnames%in%c(grep("p \\(",gsatab_colnames, value=T),grep("Stat",gsatab_colnames, value=T))])
      
      # -----------------------------------------------------------------------------
      output$loggo <- renderImage({
        filename <- "loggo.png"
        list(src = filename,
             alt = "PIANO")
      }, deleteFile = FALSE)
      
      
      # Run info --------------------------------------------------------------------
      output$text_runinfo1 <- renderUI({ 
        tmp <- gsares$info
        info <- c(
          "<h4>Gene/gene-set info:</h4>",
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
          "<h4>Analysis details:</h4>",
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
        layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,2))
        tmp <- unlist(lapply(gsares$gsc,length))
        par(mar=c(0, 4.1, 1.1, 2.1))
        h <- hist(tmp, n=50, xlim=c(0,max(tmp)), xaxt="n", xlab=NULL, main=NULL, col="forestgreen")
        par(mar=c(5, 4.1, 0, 2.1))
        boxplot(tmp, horizontal=T, ylim=c(0,max(tmp)), col="forestgreen", axes=F, xlab="Number of genes")
        axis(1, las=2)
      })
      output$gsc_boxplot2 <- renderPlot({
        layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,2))
        tmp <- unlist(lapply(gs_per_gene_list,length))
        par(mar=c(0, 4.1, 1.1, 2.1))
        h <- hist(tmp, n=50, xlim=c(0,max(tmp)), xaxt="n", xlab=NULL, main=NULL, col="mediumpurple2")
        par(mar=c(5, 4.1, 0, 2.1))
        boxplot(tmp, horizontal=T, ylim=c(0,max(tmp)), col="mediumpurple2", axes=F, xlab="Number of gene-sets")
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
        tmp <- grep("Genes",gsatab_colnames, value=T)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_padj, {
        tmp <- grep("p adj",gsatab_colnames, value=T)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_p, {
        tmp <- grep("p \\(",gsatab_colnames, value=T)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      observeEvent(input$toggle_stats, {
        tmp <- grep("Stat",gsatab_colnames, value=T)
        if(all(tmp%in%rval$visible_columns)) {
          rval$visible_columns <- rval$visible_columns[!rval$visible_columns %in% tmp]
        } else {
          rval$visible_columns <- c(rval$visible_columns,tmp)
        }
      })
      
      output$gsaresTable <- DT::renderDataTable({tmp <- rval$gs_table[gsatab_colnames[gsatab_colnames%in%rval$visible_columns]]}, 
                                                server=T, escape=F, 
                                                selection=list(mode='single', target='row'),
                                                filter="none",
                                                rownames=F,
                                                extensions=c("Scroller"),
                                                options = list(initComplete=JS(
                                                    # this is to set color format of first row
                                                    "function(settings, json) {",
                                                    "$(this.api().table().header()).css({'background-color': '#E2F3FF', 'color': '#000'});",
                                                    "}"),
                                                    dom = 'frtip',
                                                    deferRender=TRUE,
                                                    scrollY="calc(100vh - 330px)",
                                                    scrollX=T,
                                                    scroller=TRUE
                                                ))
      
      observeEvent(input$gsaresTable_rows_selected, {
        rval$sel_gs <- rval$gs_table$Name[input$gsaresTable_rows_selected]
        rval$red_gene <- ""
        updateTabsetPanel(session,"tabset1",selected="tab_gssum")
        selectRows(dataTableProxy("gsaresTable", session), list())
      })
      
      # Gene-set summary ------------------------------------------------------------
      output$text_selected_gs <- renderText({ 
        rval$sel_gs
      })
      
      observeEvent(input$link_tab_gsares, {
        updateTabsetPanel(session,"tabset1",selected="tab_gsares")
      })
      
      # Tables with info:
      output$gsTable2 <- renderTable({
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        tmp <- tmp$stats
        rownames(tmp) <- tmp[,1]
        tmp <- rbind(
        rev(c(tmp["Stat (dist.dir.dn)",2],tmp["Stat (mix.dir.dn)",2],tmp["Stat (non-dir)",2],tmp["Stat (mix.dir.up)",2],tmp["Stat (dist.dir.up)",2])),
        rev(c(tmp["p (dist.dir.dn)",2],tmp["p (mix.dir.dn)",2],tmp["p (non-dir)",2],tmp["p (mix.dir.up)",2],tmp["p (dist.dir.up)",2])),
        rev(c(tmp["p adj (dist.dir.dn)",2],tmp["p adj (mix.dir.dn)",2],tmp["p adj (non-dir)",2],tmp["p adj (mix.dir.up)",2],tmp["p adj (dist.dir.up)",2]))
        )
        rownames(tmp) <- c("Gene-set statistic","p-value","Adjusted p-value")
        colnames(tmp) <- rev(c("Distinct directional (down)","Mixed directional (down)","Non-directional","Mixed directional (up)", "Distinct directional (up)"))
        t(tmp)
      }, rownames=T)
      
      output$gsNgenes <- renderUI({
        tmp <- geneSetSummary(gsares,rval$sel_gs)
        tmp <- paste("<b>The gene-set contains ",tmp$stats[tmp$stats$Name=="Genes (tot)",2],
                     " genes </b>(",
                     tmp$stats[tmp$stats$Name=="Genes (up)",2],
                     " up, ",
                     tmp$stats[tmp$stats$Name=="Genes (dn)",2],
                     " down) <br>",
                     sep=""
        )
        HTML(tmp)
      })
      
      # Genes in gene-set info:
      output$link_view_gene <- renderUI({
        if(rval$red_gene!="") {
          list(
            actionLink("view_genes_table",paste("View genes in table")),
            HTML("<br>"),
            HTML(paste("Selected:",rval$red_gene)),
            actionLink("view_gene","(view details)"),
            HTML("<br><br>")
          )
        } else if(rval$red_gene=="") {
          list(
            actionLink("view_genes_table",paste("View genes in table")),
            HTML("<br>"),
            HTML("<em><font color='#ff0000'>You can click genes below, or in the plots to the right...</font></em><br><br>")
          )
        }
      })
      observeEvent(input$view_gene, {
        rval$sel_gene <- rval$red_gene
        updateTabsetPanel(session,"tabset1",selected="tab_geneinfo")
      }) 
      observeEvent(input$view_genes_table, {
        rval$gene_table <- genetable_all[genetable_all$`Gene ID` %in% names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats),]
        rval$sel_geneset_for_filtering <- rval$sel_gs
        updateTabsetPanel(session,"tabset1",selected="tab_genetable")
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
        sliderInput("nNeighborGenesets", NULL, 1, tmp, min(c(5,tmp)), round=T, step=1, ticks=F, width="300px")
      })
      output$gsNeighbors <- renderUI({
        tmp <- names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        tmp <- setdiff(names(gsares$gsc)[unlist(lapply(gsares$gsc, function(x) sum(tmp%in%x)>=input$nNeighborGenesets))], rval$sel_gs)
        HTML(paste(paste("<a href='#'"," onclick='Shiny.onInputChange(",'"links_genesets_click", "',tmp,'"',");'","><font color='#000000'>",tmp,"</font></a>",sep=""),collapse=", "))
      })
      observeEvent(input$filter_gs2, { # action link in ui
        tmp1 <- names(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        tmp1 <- setdiff(names(gsares$gsc)[unlist(lapply(gsares$gsc, function(x) sum(tmp1%in%x)>=input$nNeighborGenesets))], rval$sel_gs)
        tmp2 <- GSAsummaryTable(gsares)
        rval$gs_table <- tmp2[tmp2$Name %in% tmp1,]
        rval$geneset_filtering_text <- paste("Gene-sets neighboring",rval$sel_gs)
        updateTabsetPanel(session,"tabset1",selected="tab_gsares")
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
        hist(tmp2,100,col=ifelse(gsares$geneStatType=="p-signed","indianred2","grey"), add=T)
        
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
        hist(tmp2,100,col="skyblue2", add=T)
        
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
                        horizontal=T, 
                        names=c(paste("All genes down (",length(gls_all_dn),")",sep=""), 
                                paste("All genes up (",length(gls_all_up),")",sep=""),
                                "Gene-set genes (down)","Gene-set genes (up)"),
                        las=1,
                        col=c("skyblue2","indianred2","white","white"),
                        border=c("black","black","dodgerblue3","red"),
                        main="", xlab="Gene-level statistics",
                        outline=T, outcol="white")
          
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
                        horizontal=T, names=c(paste("All genes (",length(gls_all),")",sep=""),"Genes in gene-set"), las=1,
                        col=c("grey","white"),
                        main="", xlab="Gene-level statistics",
                        outline=T, outcol="white")
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
      
      output$ncharSlider <- renderUI({
        tmp <- length(geneSetSummary(gsares,rval$sel_gs)$geneLevelStats)
        sliderInput("maxNcharGenetable", NULL, 10, 1000, 50, round=T, step=1, ticks=F, width="300px")
      })
      
      output$geneTable <- DT::renderDataTable({apply(rval$gene_table,2,function(x) {
                                                                          paste("<div title='",ifelse(is.na(x),"",x),
                                                                                "' class='DT_genetable_cell'><p>",
                                                                                ifelse(is.na(x),"",substr(x,1,input$maxNcharGenetable)),
                                                                                ifelse(is.na(x),"",ifelse(nchar(x)>input$maxNcharGenetable,"...","")),
                                                                                "</div>",sep="")
                                                                        }
                                                     )
                                               }, 
        server=T, escape=F, 
        selection=list(mode='single', target='row'),
        filter="none",
        rownames=F,
        extensions=c("Scroller"),
        fillContainer=T,
        options = list(initComplete=JS(
        # this is to set color format of first row
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#E2F3FF', 'color': '#000'});",
        "}"),
        dom = 'frtip',
        deferRender=TRUE,
        scrollY="calc(100vh - 330px)",
        scrollX=T,
        scroller=TRUE)
      )
      observeEvent(input$geneTable_rows_selected, {
        rval$sel_gene <- rval$gene_table$`Gene ID`[input$geneTable_rows_selected]
        rval$red_gene <- ""
        updateTabsetPanel(session,"tabset1",selected="tab_geneinfo")
        selectRows(dataTableProxy("geneTable", session), list())
      })
      
      
      # Gene summary --------------------------------------------------------
      observeEvent(input$link_tab_genetable, { # action link in ui
        updateTabsetPanel(session,"tabset1",selected="tab_genetable")
      })
      
      output$text_selected_gene <- renderText({ 
        rval$sel_gene
      })
      
      output$info_selected_gene <- renderUI({
        HTML(paste("<b>Gene-level statistic:</b>",gsares$geneLevelStats[rval$sel_gene,],"<br>",
                   "<b>Direction:</b>",ifelse(is.na(gsares_directions[rval$sel_gene,]),"NA",ifelse(gsares_directions[rval$sel_gene,]>0,"Up","Down")),"<br>",
                   paste("<b>",colnames(genetable_all)[-c(1,2,3,4)], ":</b>",t(genetable_all[rval$sel_gene,-c(1,2,3,4)]), collapse="<br>"))
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
        updateTabsetPanel(session,"tabset1",selected="tab_gsares")
      })
      observeEvent(input$links_genesets_click, {
        rval$sel_gs <- input$links_genesets_click
        rval$red_gene <- ""
        selectRows(dataTableProxy("gsaresTable", session), list())
        updateTabsetPanel(session,"tabset1",selected="tab_gssum")
      })
      
      # Network plot --------------------------------------------------------------------
      #output$nwPlot <- renderPlot({    
      #  networkPlot(gsares, class="non")
      #})
      
      # Heatmap -------------------------------------------------------------------------
      #output$heatmap <- renderPlot({    
      #  GSAheatmap(gsares, ncharLabel=400, cex=2)
      #})
    }
  )
  
  # launch app
  if(!browser) browser <- getOption("shiny.launch.browser", interactive())
  runApp(app, launch.browser=browser)
  
}
