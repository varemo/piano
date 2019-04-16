#' Heatmap of top significant gene sets
#' 
#' This function selects the top scoring (most significant) gene sets for each
#' directionality class and produces a heatmap plot of the results.
#' 
#' This function selects the top significant gene sets in each directionality
#' class and draws a heatmap of the results. It provides a quick summary
#' alternative to the \code{\link{GSAsummaryTable}} function or the
#' \code{\link{networkPlot}}.
#' 
#' @param gsaRes an object of class \code{GSAres}, as returned from
#' \code{runGSA()}.
#' @param cutoff an integer n, so that the top n gene sets (plus possible ties)
#' in each directionality class will be included in the heatmap.
#' @param adjusted a logical, whether to use adjusted p-values or not. Note
#' that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the
#' adjusted p-values will be equal to the original p-values.
#' @param ncharLabel the number of characters to include in the row labels.
#' @param cellnote a character string selecting the information to be printed
#' inside each cell of the heatmap. Either \code{"pvalue"}, \code{"rank"},
#' \code{"nGenes"} or \code{"none"}. Note that the actual heatmap will always
#' be based on the gene set ranks.
#' @param columnnames either \code{"full"} (default) or \code{"abbr"} to use
#' full or abbreviated column labels. Will save some space for the heatmap if
#' set to \code{"abbr"}
#' @param colorkey a logical (default \code{TRUE}), whether or not to display
#' the colorkey. Will save some space for the heatmap if turned off.
#' @param colorgrad a character vector giving the color names to use in the
#' heatmap.
#' @param cex a numeric, to control the text size.
#' @return A list, returned invisibly, containing the matrix of p-values
#' (adjusted or non-adjusted depending on the settings) as represented in the
#' heatmap as well as the matrix of corresponding ranks and the matrix of
#' number of genes in each gene set (inlcuding the subset of up and down
#' regulated genes for the mixed directional classes).
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runGSA}},
#' \code{\link{GSAsummaryTable}}, \code{\link{networkPlot2}}, \code{\link{exploreGSAres}}
#' @examples
#' 
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
#'    # Make heatmap:
#'    dev.new(width=10,height=10)
#'    GSAheatmap(gsares)
#' 
GSAheatmap <- function(gsaRes, cutoff=5, adjusted=FALSE, ncharLabel=25, cellnote="pvalue", columnnames="full",
                       colorkey=TRUE, colorgrad=NULL, cex=NULL) {
   
   tmp <- try(cellnote <- match.arg(cellnote, c("pvalue","rank","nGenes","none"), several.ok=FALSE), silent=TRUE)
   if(is(tmp, "try-error")) {
      stop("argument cellnote is not valid")
   }
   if(cellnote=="pvalue") cellnote<-"medianPvalue"
   if(cellnote=="rank") cellnote<-"consensusScore"
   
   # This will allow the heatmap even with some columns NA. But need to fix this handling in 
   # consensusHeatmap and consensusScores first also....
   #if(all(is.na(gsaRes$pMixedDirUp[,1]))) gsaRes$pMixedDirUp[,1] <- rep(1,nrow(gsaRes$pMixedDirUp))
   #if(all(is.na(gsaRes$pMixedDirDn[,1]))) gsaRes$pMixedDirDn[,1] <- rep(1,nrow(gsaRes$pMixedDirDn))
   #if(all(is.na(gsaRes$pDistinctDirUp[,1]))) gsaRes$pDistinctDirUp[,1] <- rep(1,nrow(gsaRes$pDistinctDirUp))
   #if(all(is.na(gsaRes$pDistinctDirDn[,1]))) gsaRes$pDistinctDirDn[,1] <- rep(1,nrow(gsaRes$pDistinctDirDn))
   
   tmp <- FALSE
   if(all(is.na(gsaRes$pMixedDirUp[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pMixedDirDn[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pDistinctDirUp[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pDistinctDirDn[,1]))) tmp <- TRUE
   if(tmp) stop("The heatmap will only be drawn if all p-value classes exist. This is not the case for the supplied gsaRes object. Take a look at the networkPlot function instead.")
   
   res <- consensusHeatmap(resList=list(gsaRes,gsaRes), method="median", cutoff=cutoff, adjusted=adjusted, 
                     ncharLabel=ncharLabel, plot=TRUE, cellnote=cellnote, columnnames=columnnames,
                           colorkey=colorkey, colorgrad=colorgrad, cex=cex)
   invisible(res)
}
